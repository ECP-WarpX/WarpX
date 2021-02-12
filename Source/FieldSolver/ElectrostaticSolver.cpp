/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeTensorLaplacian.H>

#include <WarpX.H>

#include <memory>

using namespace amrex;

void
WarpX::ComputeSpaceChargeField (bool const reset_fields)
{
    if (reset_fields) {
        // Reset all E and B fields to 0, before calculating space-charge fields
        for (int lev = 0; lev <= max_level; lev++) {
            for (int comp=0; comp<3; comp++) {
                Efield_fp[lev][comp]->setVal(0);
                Bfield_fp[lev][comp]->setVal(0);
            }
        }
    }

    if (do_electrostatic == ElectrostaticSolverAlgo::LabFrame) {
        AddSpaceChargeFieldLabFrame();
    } else {
        // Loop over the species and add their space-charge contribution to E and B
        for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
            WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
            if (species.initialize_self_fields ||
                (do_electrostatic == ElectrostaticSolverAlgo::Relativistic)) {
                AddSpaceChargeField(species);
            }
        }
    }

}

void
WarpX::AddSpaceChargeField (WarpXParticleContainer& pc)
{

#ifdef WARPX_DIM_RZ
    amrex::Abort("The initialization of space-charge field has not yet been implemented in RZ geometry.");
#endif

    // Allocate fields for charge and potential
    const int num_levels = max_level + 1;
    Vector<std::unique_ptr<MultiFab> > rho(num_levels);
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    // Use number of guard cells used for local deposition of rho
    const int ng = guard_cells.ng_depos_rho.max();
    for (int lev = 0; lev <= max_level; lev++) {
        BoxArray nba = boxArray(lev);
        nba.surroundingNodes();
        rho[lev] = std::make_unique<MultiFab>(nba, dmap[lev], 1, ng);
        phi[lev] = std::make_unique<MultiFab>(nba, dmap[lev], 1, 1);
        phi[lev]->setVal(0.);
    }

    // Deposit particle charge density (source of Poisson solver)
    bool const local = false;
    bool const reset = true;
    bool const do_rz_volume_scaling = true;
    pc.DepositCharge(rho, local, reset, do_rz_volume_scaling);

    // Get the particle beta vector
    bool const local_average = false; // Average across all MPI ranks
    std::array<Real, 3> beta = pc.meanParticleVelocity(local_average);
    for (Real& beta_comp : beta) beta_comp /= PhysConst::c; // Normalize

    // Compute the potential phi, by solving the Poisson equation
    computePhi( rho, phi, beta, pc.self_fields_required_precision, pc.self_fields_max_iters );

    // Compute the corresponding electric and magnetic field, from the potential phi
    computeE( Efield_fp, phi, beta );
    computeB( Bfield_fp, phi, beta );

}

void
WarpX::AddSpaceChargeFieldLabFrame ()
{

#ifdef WARPX_DIM_RZ
    amrex::Abort("The calculation of space-charge field has not yet been implemented in RZ geometry.");
#endif

    // Allocate fields for charge
    // Also, zero out the phi data - is this necessary?
    const int num_levels = max_level + 1;
    Vector<std::unique_ptr<MultiFab> > rho(num_levels);
    // Use number of guard cells used for local deposition of rho
    const int ng = guard_cells.ng_depos_rho.max();
    for (int lev = 0; lev <= max_level; lev++) {
        BoxArray nba = boxArray(lev);
        nba.surroundingNodes();
        rho[lev] = std::make_unique<MultiFab>(nba, dmap[lev], 1, ng);
        rho[lev]->setVal(0.);
        phi_fp[lev]->setVal(0.);
    }

    // Deposit particle charge density (source of Poisson solver)
    for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
        WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
        bool const local = true;
        bool const reset = false;
        bool const do_rz_volume_scaling = false;
        species.DepositCharge(rho, local, reset, do_rz_volume_scaling);
    }
    for (int lev = 0; lev <= max_level; lev++) {
        ApplyFilterandSumBoundaryRho (lev, lev, *rho[lev], 0, 1);
    }

    // beta is zero in lab frame
    // Todo: use simpler finite difference form with beta=0
    std::array<Real, 3> beta = {0._rt};

    // Compute the potential phi, by solving the Poisson equation
    computePhi( rho, phi_fp, beta, self_fields_required_precision, self_fields_max_iters );

    // Compute the corresponding electric and magnetic field, from the potential phi
    computeE( Efield_fp, phi_fp, beta );
    computeB( Bfield_fp, phi_fp, beta );

}

/* Compute the potential `phi` by solving the Poisson equation with `rho` as
   a source, assuming that the source moves at a constant speed \f$\vec{\beta}\f$.
   This uses the amrex solver.

   More specifically, this solves the equation
   \f[
       \vec{\nabla}^2\phi - (\vec{\beta}\cdot\vec{\nabla})^2\phi = -\frac{\rho}{\epsilon_0}
   \f]

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computePhi (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                   amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<Real, 3> const beta,
                   Real const required_precision,
                   int const max_iters) const
{
    // Define the boundary conditions
    Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
    for (int idim=0; idim<AMREX_SPACEDIM; idim++){
        if ( Geom(0).isPeriodic(idim) ) {
            lobc[idim] = LinOpBCType::Periodic;
            hibc[idim] = LinOpBCType::Periodic;
        } else {
            // Use Dirichlet boundary condition by default.
            // Ideally, we would often want open boundary conditions here.
            lobc[idim] = LinOpBCType::Dirichlet;
            hibc[idim] = LinOpBCType::Dirichlet;

            // set the boundary potential values for the given dimension
            setDirichletBC(phi, idim);
        }
    }

    // Define the linear operator (Poisson operator)
    MLNodeTensorLaplacian linop( Geom(), boxArray(), DistributionMap() );
    linop.setDomainBC( lobc, hibc );

    // This is apparently necessary for how AMReX expects phi
    for (int lev=0; lev < rho.size(); lev++){
        phi[lev]->mult(-1.*PhysConst::ep0);
    }

    // Set the value of beta
    amrex::Array<amrex::Real,AMREX_SPACEDIM> beta_solver =
#if (AMREX_SPACEDIM==2)
        {{ beta[0], beta[2] }};  // beta_x and beta_z
#else
        {{ beta[0], beta[1], beta[2] }};
#endif
    linop.setBeta( beta_solver );

    // Solve the Poisson equation
    MLMG mlmg(linop);
    mlmg.setVerbose(2);
    mlmg.setMaxIter(max_iters);
    mlmg.solve( GetVecOfPtrs(phi), GetVecOfConstPtrs(rho), required_precision, 0.0);

    // Normalize by the correct physical constant
    for (int lev=0; lev < rho.size(); lev++){
        phi[lev]->mult(-1./PhysConst::ep0);
    }
}

/* \bried Set Dirichlet boundary conditions for the electrostatic solver.

    The given potential's values are fixed on the boundaries of the given
    dimension according to the desired values from the simulation input file,
    geometry.potential_lo and geometry.potential_hi.

   \param[inout] phi The electrostatic potential
   \param[in] idim The dimension for which the Dirichlet boundary condition is set
*/
void
WarpX::setDirichletBC(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                      const int idim) const
{
    // Get the boundary potentials specified in the simulation input file
    Vector<Real> potential_lo(AMREX_SPACEDIM);
    Vector<Real> potential_hi(AMREX_SPACEDIM);

    ParmParse pp_geom("geometry");
    pp_geom.getarr("potential_lo",potential_lo,0,AMREX_SPACEDIM);
    pp_geom.getarr("potential_hi",potential_hi,0,AMREX_SPACEDIM);

    // Print() << "Potential left = " << potential_lo[0] << "\n";
    // Print() << "Potential right = " << potential_hi[0] << "\n";

    // loop over all multigrid levels and set the boundary values
    for (int lev=0; lev <= max_level; lev++) {
        for (MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            FArrayBox& fab = (*phi[0])[mfi]; //
            Array4<Real> const& phi_gh = fab.array(); // phi include ghost cells
            const auto lo = lbound(phi_gh);
            const auto hi = ubound(phi_gh);

            //TODO use AMREX_GPU_DEVICE for acceleration
            // The boundary values are set at points lo.idim + 1 and hi.idim -1
            // since phi_gh includes ghost cells
            if (idim == 0){
                for (int k=lo.z; k<=hi.z; k++) {
                for (int j=lo.y; j<=hi.y; j++) {
                    phi_gh(lo.x+1,j,k) = potential_lo[idim];
                    phi_gh(hi.x-1,j,k) = potential_hi[idim];
                }}
            }
            if (idim == 1){
                for (int k=lo.z; k<=hi.z; k++) {
                for (int i=lo.x; i<=hi.x; i++) {
                    phi_gh(i,lo.y+1,k) = potential_lo[idim];
                    phi_gh(i,hi.y-1,k) = potential_hi[idim];
                }}
            }
            if (idim == 2){
                for (int j=lo.y; j<=hi.y; j++) {
                for (int i=lo.x; i<=hi.x; i++) {
                    phi_gh(i,j,lo.z+1) = potential_lo[idim];
                    phi_gh(i,j,hi.z-1) = potential_hi[idim];
                }}
            }
    }} // lev & MFIter

    /*
    Wei-Ting's implementation with GPU acceleration:

    // 1.0 Init V as a linear function of X
    amrex::Real Lx     = prob_hi[0] - prob_lo[0];
    amrex::Real DeltaV = potential_hi[0] - potential_lo[0];
    amrex::Real V_left = potential_lo[0];

    //### query the box size: maybe have the vars as the class var
    Vector<Real> prob_lo(AMREX_SPACEDIM);
    Vector<Real> prob_hi(AMREX_SPACEDIM);

    pp_geom.getarr("prob_lo",prob_lo,0,AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);
    pp_geom.getarr("prob_hi",prob_hi,0,AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    // does not loop into ghost cell, since phi is a nodal vector
    for (int lev=0; lev <= max_level; lev++) {
      amrex::Real dx = WarpX::CellSize(lev)[0];

      for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
         auto phi_arr   = phi[lev]->array(mfi);
         const Box& tb  = mfi.tilebox( phi[lev]->ixType().toIntVect() );

         amrex::ParallelFor( tb,
             [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                phi_arr(i,j,k) = (i*dx/Lx)*DeltaV + V_left;
             } // loop ijk
         );
    }} // lev & MRIter
    */

}

/* \bried Compute the electric field that corresponds to `phi`, and
          add it to the set of MultiFab `E`.

   The electric field is calculated by assuming that the source that
   produces the `phi` potential is moving with a constant speed \f$\vec{\beta}\f$:
   \f[
    \vec{E} = -\vec{\nabla}\phi + (\vec{\beta}\cdot\vec{\beta})\phi \vec{\beta}
   \f]
   (where the second term represent the term \f$\partial_t \vec{A}\f$, in
    the case of a moving source)

   \param[inout] E Electric field on the grid
   \param[in] phi The potential from which to compute the electric field
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computeE (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                 const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                 std::array<amrex::Real, 3> const beta ) const
{
    for (int lev = 0; lev <= max_level; lev++) {

        const Real* dx = Geom(lev).CellSize();

#ifdef _OPENMP
#    pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Real inv_dx = 1./dx[0];
#if (AMREX_SPACEDIM == 3)
            const Real inv_dy = 1./dx[1];
            const Real inv_dz = 1./dx[2];
#else
            const Real inv_dz = 1./dx[1];
#endif
            const Box& tbx  = mfi.tilebox( E[lev][0]->ixType().toIntVect() );
#if (AMREX_SPACEDIM == 3)
            const Box& tby  = mfi.tilebox( E[lev][1]->ixType().toIntVect() );
#endif
            const Box& tbz  = mfi.tilebox( E[lev][2]->ixType().toIntVect() );

            const auto& phi_arr = phi[lev]->array(mfi);
            const auto& Ex_arr = (*E[lev][0])[mfi].array();
#if (AMREX_SPACEDIM == 3)
            const auto& Ey_arr = (*E[lev][1])[mfi].array();
#endif
            const auto& Ez_arr = (*E[lev][2])[mfi].array();

            Real beta_x = beta[0];
            Real beta_y = beta[1];
            Real beta_z = beta[2];

            // Calculate the electric field
            // Use discretized derivative that matches the staggering of the grid.
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex_arr(i,j,k) +=
                        +(beta_x*beta_x-1)*inv_dx*( phi_arr(i+1,j,k)-phi_arr(i,j,k) )
                        +beta_x*beta_y*0.25*inv_dy*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j-1,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i+1,j-1,k))
                        +beta_x*beta_z*0.25*inv_dz*(phi_arr(i  ,j,k+1)-phi_arr(i  ,j,k-1)
                                                  + phi_arr(i+1,j,k+1)-phi_arr(i+1,j,k-1));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ey_arr(i,j,k) +=
                        +beta_y*beta_x*0.25*inv_dx*(phi_arr(i+1,j  ,k)-phi_arr(i-1,j  ,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i-1,j+1,k))
                        +(beta_y*beta_y-1)*inv_dy*( phi_arr(i,j+1,k)-phi_arr(i,j,k) )
                        +beta_y*beta_z*0.25*inv_dz*(phi_arr(i,j  ,k+1)-phi_arr(i,j  ,k-1)
                                                  + phi_arr(i,j+1,k+1)-phi_arr(i,j+1,k-1));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) +=
                        +beta_z*beta_x*0.25*inv_dx*(phi_arr(i+1,j,k  )-phi_arr(i-1,j,k  )
                                                  + phi_arr(i+1,j,k+1)-phi_arr(i-1,j,k+1))
                        +beta_z*beta_y*0.25*inv_dy*(phi_arr(i,j+1,k  )-phi_arr(i,j-1,k  )
                                                  + phi_arr(i,j+1,k+1)-phi_arr(i,j-1,k+1))
                        +(beta_y*beta_z-1)*inv_dz*( phi_arr(i,j,k+1)-phi_arr(i,j,k) );
                }
            );
#else
            amrex::ParallelFor( tbx, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex_arr(i,j,k) +=
                        +(beta_x*beta_x-1)*inv_dx*( phi_arr(i+1,j,k)-phi_arr(i,j,k) )
                        +beta_x*beta_z*0.25*inv_dz*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j-1,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i+1,j-1,k));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) +=
                        +beta_z*beta_x*0.25*inv_dx*(phi_arr(i+1,j  ,k)-phi_arr(i-1,j  ,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i-1,j+1,k))
                        +(beta_y*beta_z-1)*inv_dz*( phi_arr(i,j+1,k)-phi_arr(i,j,k) );
                }
            );
#endif
        }
    }
}


/* \bried Compute the magnetic field that corresponds to `phi`, and
          add it to the set of MultiFab `B`.

   The magnetic field is calculated by assuming that the source that
   produces the `phi` potential is moving with a constant speed \f$\vec{\beta}\f$:
   \f[
    \vec{B} = -\frac{1}{c}\vec{\beta}\times\vec{\nabla}\phi
   \f]
   (this represents the term \f$\vec{\nabla} \times \vec{A}\f$, in the case of a moving source)

   \param[inout] E Electric field on the grid
   \param[in] phi The potential from which to compute the electric field
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computeB (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& B,
                 const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                 std::array<amrex::Real, 3> const beta ) const
{
    for (int lev = 0; lev <= max_level; lev++) {

        const Real* dx = Geom(lev).CellSize();

#ifdef _OPENMP
#    pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Real inv_dx = 1./dx[0];
#if (AMREX_SPACEDIM == 3)
            const Real inv_dy = 1./dx[1];
            const Real inv_dz = 1./dx[2];
#else
            const Real inv_dz = 1./dx[1];
#endif
            const Box& tbx  = mfi.tilebox( B[lev][0]->ixType().toIntVect() );
            const Box& tby  = mfi.tilebox( B[lev][1]->ixType().toIntVect() );
            const Box& tbz  = mfi.tilebox( B[lev][2]->ixType().toIntVect() );

            const auto& phi_arr = phi[0]->array(mfi);
            const auto& Bx_arr = (*B[lev][0])[mfi].array();
            const auto& By_arr = (*B[lev][1])[mfi].array();
            const auto& Bz_arr = (*B[lev][2])[mfi].array();

            Real beta_x = beta[0];
            Real beta_y = beta[1];
            Real beta_z = beta[2];

            constexpr Real inv_c = 1./PhysConst::c;

            // Calculate the magnetic field
            // Use discretized derivative that matches the staggering of the grid.
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bx_arr(i,j,k) += inv_c * (
                        -beta_y*inv_dz*0.5*(phi_arr(i,j  ,k+1)-phi_arr(i,j  ,k)
                                          + phi_arr(i,j+1,k+1)-phi_arr(i,j+1,k))
                        +beta_z*inv_dy*0.5*(phi_arr(i,j+1,k  )-phi_arr(i,j,k  )
                                          + phi_arr(i,j+1,k+1)-phi_arr(i,j,k+1)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    By_arr(i,j,k) += inv_c * (
                        -beta_z*inv_dx*0.5*(phi_arr(i+1,j,k  )-phi_arr(i,j,k  )
                                          + phi_arr(i+1,j,k+1)-phi_arr(i,j,k+1))
                        +beta_x*inv_dz*0.5*(phi_arr(i  ,j,k+1)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j,k+1)-phi_arr(i+1,j,k)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bz_arr(i,j,k) += inv_c * (
                        -beta_x*inv_dy*0.5*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i+1,j,k))
                        +beta_y*inv_dx*0.5*(phi_arr(i+1,j  ,k)-phi_arr(i,j  ,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i,j+1,k)));
                }
            );
#else
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bx_arr(i,j,k) += inv_c * (
                        -beta_y*inv_dz*( phi_arr(i,j+1,k)-phi_arr(i,j,k) ));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    By_arr(i,j,k) += inv_c * (
                        -beta_z*inv_dx*0.5*(phi_arr(i+1,j  ,k)-phi_arr(i,j  ,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i,j+1,k))
                        +beta_x*inv_dz*0.5*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i+1,j,k)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bz_arr(i,j,k) += inv_c * (
                        +beta_y*inv_dx*( phi_arr(i+1,j,k)-phi_arr(i,j,k) ));
                }
            );
#endif
        }
    }
}
