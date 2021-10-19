/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "FieldSolver/ElectrostaticSolver.H"
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Python/WarpX_py.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MFIter.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <array>
#include <memory>
#include <string>

using namespace amrex;

void
WarpX::ComputeSpaceChargeField (bool const reset_fields)
{
    WARPX_PROFILE("WarpX::ComputeSpaceChargeField");
    if (reset_fields) {
        // Reset all E and B fields to 0, before calculating space-charge fields
        WARPX_PROFILE("WarpX::ComputeSpaceChargeField::reset_fields");
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
    // Transfer fields from 'fp' array to 'aux' array.
    // This is needed when using momentum conservation
    // since they are different arrays in that case.
    UpdateAuxilaryData();
    FillBoundaryAux(guard_cells.ng_UpdateAux);

}

void
WarpX::AddSpaceChargeField (WarpXParticleContainer& pc)
{
    WARPX_PROFILE("WarpX::AddSpaceChargeField");

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!field_boundary_handler.bcs_set) field_boundary_handler.definePhiBCs();

#ifdef WARPX_DIM_RZ
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_rz_azimuthal_modes == 1,
                                     "Error: RZ electrostatic only implemented for a single mode");
#endif

    // Allocate fields for charge and potential
    const int num_levels = max_level + 1;
    Vector<std::unique_ptr<MultiFab> > rho(num_levels);
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    // Use number of guard cells used for local deposition of rho
    const amrex::IntVect ng = guard_cells.ng_depos_rho;
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
    computePhi( rho, phi, beta, pc.self_fields_required_precision, pc.self_fields_max_iters, pc.self_fields_verbosity );

    // Compute the corresponding electric and magnetic field, from the potential phi
    computeE( Efield_fp, phi, beta );
    computeB( Bfield_fp, phi, beta );

}

void
WarpX::AddSpaceChargeFieldLabFrame ()
{
    WARPX_PROFILE("WarpX::AddSpaceChargeFieldLabFrame");

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!field_boundary_handler.bcs_set) field_boundary_handler.definePhiBCs();

#ifdef WARPX_DIM_RZ
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_rz_azimuthal_modes == 1,
                                     "Error: RZ electrostatic only implemented for a single mode");
#endif

    // reset rho_fp before depositing charge density for this step
    for (int lev = 0; lev <= max_level; lev++) {
        rho_fp[lev]->setVal(0.);
    }

    // Deposit particle charge density (source of Poisson solver)
    bool const local = true;
    bool const interpolate_across_levels = false;
    bool const reset = false;
    bool const do_rz_volume_scaling = false;
    for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
        WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
        species.DepositCharge(
            rho_fp, local, reset, do_rz_volume_scaling, interpolate_across_levels
        );
    }
#ifdef WARPX_DIM_RZ
    for (int lev = 0; lev <= max_level; lev++) {
        ApplyInverseVolumeScalingToChargeDensity(rho_fp[lev].get(), lev);
    }
#endif
    SyncRho(); // Apply filter, perform MPI exchange, interpolate across levels

    // beta is zero in lab frame
    // Todo: use simpler finite difference form with beta=0
    std::array<Real, 3> beta = {0._rt};

    // Compute the potential phi, by solving the Poisson equation
    if (warpx_py_poissonsolver) warpx_py_poissonsolver();
    else computePhi( rho_fp, phi_fp, beta, self_fields_required_precision,
                     self_fields_max_iters, self_fields_verbosity );

    // Compute the electric field. Note that if an EB is used the electric
    // field will be calculated in the computePhi call.
#ifndef AMREX_USE_EB
    computeE( Efield_fp, phi_fp, beta );
#else
    if (warpx_py_poissonsolver) computeE( Efield_fp, phi_fp, beta );
#endif

    // Compute the magnetic field
    computeB( Bfield_fp, phi_fp, beta );
}

/* Compute the potential `phi` by solving the Poisson equation with `rho` as
   a source, assuming that the source moves at a constant speed \f$\vec{\beta}\f$.
   This uses the amrex solver.

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computePhi (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                   amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<Real, 3> const beta,
                   Real const required_precision,
                   int const max_iters,
                   int const verbosity) const
{
#ifdef WARPX_DIM_RZ
    computePhiRZ( rho, phi, beta, required_precision, max_iters, verbosity );
#else
    computePhiCartesian( rho, phi, beta, required_precision, max_iters, verbosity );
#endif

}

#ifdef WARPX_DIM_RZ
/* Compute the potential `phi` in cylindrical geometry by solving the Poisson equation
   with `rho` as a source, assuming that the source moves at a constant
   speed \f$\vec{\beta}\f$.
   This uses the amrex solver.

   More specifically, this solves the equation
   \f[
       \vec{\nabla}^2 r \phi - (\vec{\beta}\cdot\vec{\nabla})^2 r \phi = -\frac{r \rho}{\epsilon_0}
   \f]

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
   \param[in] beta Represents the velocity of the source of `phi`
*/
void
WarpX::computePhiRZ (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                   amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<Real, 3> const beta,
                   Real const required_precision,
                   int const max_iters,
                   int const verbosity) const
{
    // Create a new geometry with the z coordinate scaled by gamma
    amrex::Real const gamma = std::sqrt(1._rt/(1._rt - beta[2]*beta[2]));

    amrex::Vector<amrex::Geometry> geom_scaled(max_level + 1);
    for (int lev = 0; lev <= max_level; ++lev) {
        const amrex::Geometry & geom_lev = Geom(lev);
        const amrex::Real* current_lo = geom_lev.ProbLo();
        const amrex::Real* current_hi = geom_lev.ProbHi();
        amrex::Real scaled_lo[AMREX_SPACEDIM];
        amrex::Real scaled_hi[AMREX_SPACEDIM];
        scaled_lo[0] = current_lo[0];
        scaled_hi[0] = current_hi[0];
        scaled_lo[1] = current_lo[1]*gamma;
        scaled_hi[1] = current_hi[1]*gamma;
        amrex::RealBox rb = RealBox(scaled_lo, scaled_hi);
        geom_scaled[lev].define(geom_lev.Domain(), &rb);
    }

    // Setup the sigma = radius
    // sigma must be cell centered
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > sigma(max_level+1);
    for (int lev = 0; lev <= max_level; ++lev) {
        const amrex::Real * problo = geom_scaled[lev].ProbLo();
        const amrex::Real * dx = geom_scaled[lev].CellSize();
        const amrex::Real rmin = problo[0];
        const amrex::Real dr = dx[0];

        amrex::BoxArray nba = boxArray(lev);
        nba.enclosedCells(); // Get cell centered array (correct?)
        sigma[lev] = std::make_unique<MultiFab>(nba, dmap[lev], 1, 0);
        for ( MFIter mfi(*sigma[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const amrex::Box& tbx = mfi.tilebox();
            const amrex::Dim3 lo = amrex::lbound(tbx);
            const int irmin = lo.x;
            Array4<amrex::Real> const& sigma_arr = sigma[lev]->array(mfi);
            amrex::ParallelFor( tbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) {
                    sigma_arr(i,j,0) = rmin + (i - irmin + 0.5_rt)*dr;
                }
            );
        }

        // Also, multiply rho by radius (rho is node centered)
        // Note that this multiplication is not undone since rho is
        // a temporary array.
        for ( MFIter mfi(*rho[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const amrex::Box& tbx = mfi.tilebox();
            const amrex::Dim3 lo = amrex::lbound(tbx);
            const int irmin = lo.x;
            int const ncomp = rho[lev]->nComp(); // This should be 1!
            Array4<Real> const& rho_arr = rho[lev]->array(mfi);
            amrex::ParallelFor(tbx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/, int icomp)
            {
                amrex::Real r = rmin + (i - irmin)*dr;
                if (r == 0.) {
                    // dr/3 is used to be consistent with the finite volume formulism
                    // that is used to solve Poisson's equation
                    rho_arr(i,j,0,icomp) *= dr/3._rt;
                } else {
                    rho_arr(i,j,0,icomp) *= r;
                }
            });
        }
    }

    // get the potential at the current time
    amrex::Array<amrex::Real,AMREX_SPACEDIM> phi_bc_values_lo;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> phi_bc_values_hi;
    phi_bc_values_lo[1] = field_boundary_handler.potential_zlo(gett_new(0));
    phi_bc_values_hi[1] = field_boundary_handler.potential_zhi(gett_new(0));

    // set the boundary potential values if needed
    setPhiBC(phi, phi_bc_values_lo, phi_bc_values_hi);

    // Define the linear operator (Poisson operator)
    MLNodeLaplacian linop( geom_scaled, boxArray(), dmap );

    for (int lev = 0; lev <= max_level; ++lev) {
        linop.setSigma( lev, *sigma[lev] );
    }

    for (int lev=0; lev < rho.size(); lev++){
        rho[lev]->mult(-1._rt/PhysConst::ep0);
    }

    // Solve the Poisson equation
    linop.setDomainBC( field_boundary_handler.lobc, field_boundary_handler.hibc );
    MLMG mlmg(linop);
    mlmg.setVerbose(verbosity);
    mlmg.setMaxIter(max_iters);
    mlmg.solve( GetVecOfPtrs(phi), GetVecOfConstPtrs(rho), required_precision, 0.0);
}

#else

/* Compute the potential `phi` in Cartesian geometry by solving the Poisson equation
   with `rho` as a source, assuming that the source moves at a constant
   speed \f$\vec{\beta}\f$.
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
WarpX::computePhiCartesian (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                            std::array<Real, 3> const beta,
                            Real const required_precision,
                            int const max_iters,
                            int const verbosity) const
{
    // get the potential at the current time
    amrex::Array<amrex::Real,AMREX_SPACEDIM> phi_bc_values_lo;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> phi_bc_values_hi;
    phi_bc_values_lo[0] = field_boundary_handler.potential_xlo(gett_new(0));
    phi_bc_values_hi[0] = field_boundary_handler.potential_xhi(gett_new(0));
#if (AMREX_SPACEDIM==2)
    phi_bc_values_lo[1] = field_boundary_handler.potential_zlo(gett_new(0));
    phi_bc_values_hi[1] = field_boundary_handler.potential_zhi(gett_new(0));
#elif (AMREX_SPACEDIM==3)
    phi_bc_values_lo[1] = field_boundary_handler.potential_ylo(gett_new(0));
    phi_bc_values_hi[1] = field_boundary_handler.potential_yhi(gett_new(0));
    phi_bc_values_lo[2] = field_boundary_handler.potential_zlo(gett_new(0));
    phi_bc_values_hi[2] = field_boundary_handler.potential_zhi(gett_new(0));
#endif

    setPhiBC(phi, phi_bc_values_lo, phi_bc_values_hi);

#ifndef AMREX_USE_EB
    // Define the linear operator (Poisson operator)
    MLNodeTensorLaplacian linop( Geom(), boxArray(), DistributionMap() );

    // Set the value of beta
    amrex::Array<amrex::Real,AMREX_SPACEDIM> beta_solver =
#   if (AMREX_SPACEDIM==2)
        {{ beta[0], beta[2] }};  // beta_x and beta_z
#   else
        {{ beta[0], beta[1], beta[2] }};
#   endif
    linop.setBeta( beta_solver );

#else

    // With embedded boundary: extract EB info
    LPInfo info;
    Vector<EBFArrayBoxFactory const*> eb_factory;
    eb_factory.resize(max_level+1);
    for (int lev = 0; lev <= max_level; ++lev) {
      eb_factory[lev] = &WarpX::fieldEBFactory(lev);
    }
    MLEBNodeFDLaplacian linop( Geom(), boxArray(), dmap, info, eb_factory);

    // Note: this assumes that the beam is propagating along
    // one of the axes of the grid, i.e. that only *one* of the Cartesian
    // components of `beta` is non-negligible.
    linop.setSigma({AMREX_D_DECL(
        1._rt-beta[0]*beta[0], 1._rt-beta[1]*beta[1], 1._rt-beta[2]*beta[2])});

    // if the EB potential only depends on time, the potential can be passed
    // as a float instead of a callable
    if (field_boundary_handler.phi_EB_only_t) {
        linop.setEBDirichlet(field_boundary_handler.potential_eb_t(gett_new(0)));
    }
    else linop.setEBDirichlet(field_boundary_handler.getPhiEB(gett_new(0)));
#endif

    // Solve the Poisson equation
    linop.setDomainBC( field_boundary_handler.lobc, field_boundary_handler.hibc );

    for (int lev=0; lev < rho.size(); lev++){
        rho[lev]->mult(-1._rt/PhysConst::ep0);
    }

    MLMG mlmg(linop);
    mlmg.setVerbose(verbosity);
    mlmg.setMaxIter(max_iters);
    mlmg.solve( GetVecOfPtrs(phi), GetVecOfConstPtrs(rho), required_precision, 0.0);

#ifdef AMREX_USE_EB
    // use amrex to directly calculate the electric field since with EB's the
    // simple finite difference scheme in WarpX::computeE sometimes fails
    if (do_electrostatic == ElectrostaticSolverAlgo::LabFrame)
    {
        for (int lev = 0; lev <= max_level; ++lev) {
#if (AMREX_SPACEDIM==2)
            mlmg.getGradSolution(
                {amrex::Array<amrex::MultiFab*,2>{
                    get_pointer_Efield_fp(lev, 0),get_pointer_Efield_fp(lev, 2)
                    }}
            );
#elif (AMREX_SPACEDIM==3)
            mlmg.getGradSolution(
                {amrex::Array<amrex::MultiFab*,3>{
                    get_pointer_Efield_fp(lev, 0),get_pointer_Efield_fp(lev, 1),
                    get_pointer_Efield_fp(lev, 2)
                    }}
            );
            get_pointer_Efield_fp(lev, 1)->mult(-1._rt);
#endif
            get_pointer_Efield_fp(lev, 0)->mult(-1._rt);
            get_pointer_Efield_fp(lev, 2)->mult(-1._rt);
        }
    }
#endif
}
#endif

/* \brief Set Dirichlet boundary conditions for the electrostatic solver.

    The given potential's values are fixed on the boundaries of the given
    dimension according to the desired values from the simulation input file,
    boundary.potential_lo and boundary.potential_hi.

   \param[inout] phi The electrostatic potential
   \param[in] idim The dimension for which the Dirichlet boundary condition is set
*/
void
WarpX::setPhiBC( amrex::Vector<std::unique_ptr<amrex::MultiFab>>& phi,
                 Array<amrex::Real,AMREX_SPACEDIM>& phi_bc_values_lo,
                 Array<amrex::Real,AMREX_SPACEDIM>& phi_bc_values_hi ) const
{
    // check if any dimension has non-periodic boundary conditions
    if (!field_boundary_handler.has_non_periodic) return;

    auto dirichlet_flag = field_boundary_handler.dirichlet_flag;

    // loop over all mesh refinement levels and set the boundary values
    for (int lev=0; lev <= max_level; lev++) {

        amrex::Box domain = Geom(lev).Domain();
        domain.surroundingNodes();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            // Extract the potential
            auto phi_arr = phi[lev]->array(mfi);
            // Extract tileboxes for which to loop
            const Box& tb  = mfi.tilebox( phi[lev]->ixType().toIntVect() );

            // loop over dimensions
            for (int idim=0; idim<AMREX_SPACEDIM; idim++){
                // check if neither boundaries in this dimension should be set
                if (!(dirichlet_flag[2*idim] || dirichlet_flag[2*idim+1])) continue;

                // a check can be added below to test if the boundary values
                // are already correct, in which case the ParallelFor over the
                // cells can be skipped

                if (!domain.strictly_contains(tb)) {
                    amrex::ParallelFor( tb,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                            IntVect iv(AMREX_D_DECL(i,j,k));

                            if (dirichlet_flag[2*idim] && iv[idim] == domain.smallEnd(idim)){
                                phi_arr(i,j,k) = phi_bc_values_lo[idim];
                            }
                            if (dirichlet_flag[2*idim+1] && iv[idim] == domain.bigEnd(idim)) {
                                phi_arr(i,j,k) = phi_bc_values_hi[idim];
                            }

                        } // loop ijk
                    );
                }
            } // idim
    }} // lev & MFIter
}

/* \brief Compute the electric field that corresponds to `phi`, and
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

#ifdef AMREX_USE_OMP
#    pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Real inv_dx = 1._rt/dx[0];
#if (AMREX_SPACEDIM == 3)
            const Real inv_dy = 1._rt/dx[1];
            const Real inv_dz = 1._rt/dx[2];
#else
            const Real inv_dz = 1._rt/dx[1];
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
                        +beta_x*beta_y*0.25_rt*inv_dy*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j-1,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i+1,j-1,k))
                        +beta_x*beta_z*0.25_rt*inv_dz*(phi_arr(i  ,j,k+1)-phi_arr(i  ,j,k-1)
                                                  + phi_arr(i+1,j,k+1)-phi_arr(i+1,j,k-1));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ey_arr(i,j,k) +=
                        +beta_y*beta_x*0.25_rt*inv_dx*(phi_arr(i+1,j  ,k)-phi_arr(i-1,j  ,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i-1,j+1,k))
                        +(beta_y*beta_y-1)*inv_dy*( phi_arr(i,j+1,k)-phi_arr(i,j,k) )
                        +beta_y*beta_z*0.25_rt*inv_dz*(phi_arr(i,j  ,k+1)-phi_arr(i,j  ,k-1)
                                                  + phi_arr(i,j+1,k+1)-phi_arr(i,j+1,k-1));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) +=
                        +beta_z*beta_x*0.25_rt*inv_dx*(phi_arr(i+1,j,k  )-phi_arr(i-1,j,k  )
                                                  + phi_arr(i+1,j,k+1)-phi_arr(i-1,j,k+1))
                        +beta_z*beta_y*0.25_rt*inv_dy*(phi_arr(i,j+1,k  )-phi_arr(i,j-1,k  )
                                                  + phi_arr(i,j+1,k+1)-phi_arr(i,j-1,k+1))
                        +(beta_y*beta_z-1)*inv_dz*( phi_arr(i,j,k+1)-phi_arr(i,j,k) );
                }
            );
#else
            amrex::ParallelFor( tbx, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex_arr(i,j,k) +=
                        +(beta_x*beta_x-1)*inv_dx*( phi_arr(i+1,j,k)-phi_arr(i,j,k) )
                        +beta_x*beta_z*0.25_rt*inv_dz*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j-1,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i+1,j-1,k));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez_arr(i,j,k) +=
                        +beta_z*beta_x*0.25_rt*inv_dx*(phi_arr(i+1,j  ,k)-phi_arr(i-1,j  ,k)
                                                  + phi_arr(i+1,j+1,k)-phi_arr(i-1,j+1,k))
                        +(beta_y*beta_z-1)*inv_dz*( phi_arr(i,j+1,k)-phi_arr(i,j,k) );
                }
            );
#endif
        }
    }
}


/* \brief Compute the magnetic field that corresponds to `phi`, and
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

#ifdef AMREX_USE_OMP
#    pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Real inv_dx = 1._rt/dx[0];
#if (AMREX_SPACEDIM == 3)
            const Real inv_dy = 1._rt/dx[1];
            const Real inv_dz = 1._rt/dx[2];
#else
            const Real inv_dz = 1._rt/dx[1];
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

            constexpr Real inv_c = 1._rt/PhysConst::c;

            // Calculate the magnetic field
            // Use discretized derivative that matches the staggering of the grid.
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor( tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bx_arr(i,j,k) += inv_c * (
                        -beta_y*inv_dz*0.5_rt*(phi_arr(i,j  ,k+1)-phi_arr(i,j  ,k)
                                          + phi_arr(i,j+1,k+1)-phi_arr(i,j+1,k))
                        +beta_z*inv_dy*0.5_rt*(phi_arr(i,j+1,k  )-phi_arr(i,j,k  )
                                          + phi_arr(i,j+1,k+1)-phi_arr(i,j,k+1)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    By_arr(i,j,k) += inv_c * (
                        -beta_z*inv_dx*0.5_rt*(phi_arr(i+1,j,k  )-phi_arr(i,j,k  )
                                          + phi_arr(i+1,j,k+1)-phi_arr(i,j,k+1))
                        +beta_x*inv_dz*0.5_rt*(phi_arr(i  ,j,k+1)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j,k+1)-phi_arr(i+1,j,k)));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Bz_arr(i,j,k) += inv_c * (
                        -beta_x*inv_dy*0.5_rt*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i+1,j,k))
                        +beta_y*inv_dx*0.5_rt*(phi_arr(i+1,j  ,k)-phi_arr(i,j  ,k)
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
                        -beta_z*inv_dx*0.5_rt*(phi_arr(i+1,j  ,k)-phi_arr(i,j  ,k)
                                          + phi_arr(i+1,j+1,k)-phi_arr(i,j+1,k))
                        +beta_x*inv_dz*0.5_rt*(phi_arr(i  ,j+1,k)-phi_arr(i  ,j,k)
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

void ElectrostaticSolver::BoundaryHandler::definePhiBCs ( )
{
#ifdef WARPX_DIM_RZ
    lobc[0] = LinOpBCType::Neumann;
    hibc[0] = LinOpBCType::Dirichlet;
    dirichlet_flag[0] = false;
    dirichlet_flag[1] = false;
    int dim_start=1;
#else
    int dim_start=0;
#endif
    for (int idim=dim_start; idim<AMREX_SPACEDIM; idim++){
        if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::Periodic
             && WarpX::field_boundary_hi[idim] == FieldBoundaryType::Periodic ) {
            lobc[idim] = LinOpBCType::Periodic;
            hibc[idim] = LinOpBCType::Periodic;
            dirichlet_flag[idim*2] = false;
            dirichlet_flag[idim*2+1] = false;
        }
        else {
            has_non_periodic = true;
            if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC ) {
                lobc[idim] = LinOpBCType::Dirichlet;
                dirichlet_flag[idim*2] = true;
            }
            else if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::None ) {
                lobc[idim] = LinOpBCType::Neumann;
                dirichlet_flag[idim*2] = false;
            }
            else {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                    "Field boundary conditions have to be either periodic, PEC or none "
                    "when using the electrostatic solver"
                );
            }

            if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC ) {
                hibc[idim] = LinOpBCType::Dirichlet;
                dirichlet_flag[idim*2+1] = true;
            }
            else if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::None ) {
                hibc[idim] = LinOpBCType::Neumann;
                dirichlet_flag[idim*2+1] = false;
            }
            else {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                    "Field boundary conditions have to be either periodic, PEC or none "
                    "when using the electrostatic solver"
                );
            }
        }
    }
    bcs_set = true;
}

void ElectrostaticSolver::BoundaryHandler::buildParsers ()
{
    potential_xlo_parser = makeParser(potential_xlo_str, {"t"});
    potential_xhi_parser = makeParser(potential_xhi_str, {"t"});
    potential_ylo_parser = makeParser(potential_ylo_str, {"t"});
    potential_yhi_parser = makeParser(potential_yhi_str, {"t"});
    potential_zlo_parser = makeParser(potential_zlo_str, {"t"});
    potential_zhi_parser = makeParser(potential_zhi_str, {"t"});
    potential_eb_parser = makeParser(potential_eb_str, {"x", "y", "z", "t"});

    potential_xlo = potential_xlo_parser.compile<1>();
    potential_xhi = potential_xhi_parser.compile<1>();
    potential_ylo = potential_ylo_parser.compile<1>();
    potential_yhi = potential_yhi_parser.compile<1>();
    potential_zlo = potential_zlo_parser.compile<1>();
    potential_zhi = potential_zhi_parser.compile<1>();

    // check if the EB potential is a function of space or only of time
    std::set<std::string> eb_symbols = potential_eb_parser.symbols();
    if ((eb_symbols.count("x") != 0) || (eb_symbols.count("y") != 0)
            || (eb_symbols.count("z") != 0)) {
        potential_eb = potential_eb_parser.compile<4>();
        phi_EB_only_t = false;
    }
    else {
        potential_eb_parser = makeParser(potential_eb_str, {"t"});
        potential_eb_t = potential_eb_parser.compile<1>();
    }
}
