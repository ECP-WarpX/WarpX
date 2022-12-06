/* Copyright 2022 S. Eric Clark, LLNL
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "FieldSolver/MagnetostaticSolver/MagnetostaticSolver.H"
#include "FieldSolver/MagnetostaticSolver/MagnetostaticSolver_K.H"
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Python/WarpX_py.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Parallelization/WarpXComm_K.H"

#include <ablastr/utils/Communication.H>
#include <ablastr/warn_manager/WarnManager.H>

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
#include <AMReX_MFInterp_C.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EBFabFactory.H>
#endif

#include <array>
#include <memory>
#include <string>

using namespace amrex;

void
WarpX::ComputeCurrentDensityField()
{
    WARPX_PROFILE("WarpX::ComputeCurrentDensityField");
    // Fields have been reset in Electrostatic solver for this time step, these fields
    // are added into the E & B fields after electrostatic solve

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(do_current_centering == true,
                                     "Error: Magnetostatic solver requires warpx.do_current_centering=true");

    AddCurrentDensityFieldLabFrame();

    // Transfer fields from 'fp' array to 'aux' array.
    // This is needed when using momentum conservation
    // since they are different arrays in that case.
    UpdateAuxilaryData();
    FillBoundaryAux(guard_cells.ng_UpdateAux);

}

void
WarpX::AddCurrentDensityFieldLabFrame()
{
    WARPX_PROFILE("WarpX::AddCurrentDensityFieldLabFrame");

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!m_vector_poisson_boundary_handler.bcs_set) {
        m_vector_poisson_boundary_handler.defineVectorPotentialBCs();
    }

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(n_rz_azimuthal_modes == 1,
                                     "Error: RZ magnetostatic only implemented for a single mode");
#endif

    // reset current_fp before depositing current density for this step
    for (int lev = 0; lev <= max_level; lev++) {
        for (int dim=0; dim < 3; dim++) {
            current_fp_nodal[lev][dim]->setVal(0.);
        }
    }

    // Deposit current density (source of Poisson solver)
    for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
        WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
        if (!species.do_not_deposit) {
            species.DepositCurrent(current_fp_nodal, dt[0], 0.);
        }
    }

    // This will interpolate to staggered grid, filter, and synchronize
    SyncCurrent(current_fp, current_cp); // Apply filter, perform MPI exchange, interpolate across levels

#ifdef WARPX_DIM_RZ
    for (int lev = 0; lev <= max_level; lev++) {
        ApplyInverseVolumeScalingToCurrentDensity(current_fp[lev][0].get(),
                                                  current_fp[lev][1].get(),
                                                  current_fp[lev][2].get(), lev);
    }
#endif

    // NOTE:  (SEC 12/2/22)
    // Interpolate fields back to nodal current
    // This could be done faster, but want to check if answer is right
    WarpX &warpx = WarpX::GetInstance();

    // Grab Interpolation Coefficients
    // Order of finite-order centering of fields
    const int fg_nox = WarpX::field_centering_nox;
    const int fg_noy = WarpX::field_centering_noy;
    const int fg_noz = WarpX::field_centering_noz;

    // Device vectors of stencil coefficients used for finite-order centering of fields
    amrex::Real const * stencil_coeffs_x = warpx.device_field_centering_stencil_coeffs_x.data();
    amrex::Real const * stencil_coeffs_y = warpx.device_field_centering_stencil_coeffs_y.data();
    amrex::Real const * stencil_coeffs_z = warpx.device_field_centering_stencil_coeffs_z.data();

    for (int lev=0; lev<=finest_level; lev++)
    {
        for (int idim = 0; idim<3; idim++)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*(current_fp_nodal[lev][idim]), TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                IntVect const src_stag = current_fp[lev][idim]->ixType().toIntVect();
                IntVect const dst_stag = current_fp_nodal[lev][idim]->ixType().toIntVect();

                Array4<amrex::Real const> const& src_arr = current_fp[lev][idim]->const_array(mfi);
                Array4<amrex::Real> const& dst_arr = current_fp_nodal[lev][idim]->array(mfi);

                // Loop includes ghost cells (`growntilebox`)
                // (input arrays will be padded with zeros beyond ghost cells
                // for out-of-bound accesses due to large-stencil operations)
                Box bx = mfi.tilebox(dst_stag, current_fp_nodal[lev][idim]->nGrowVect());

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, dst_arr, src_arr, dst_stag, src_stag, fg_nox, fg_noy, fg_noz,
                        stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
                });
            }

            const amrex::Periodicity& curr_period = geom[lev].periodicity();
            // Synchronize the ghost cells, do halo exchange
            ablastr::utils::communication::FillBoundary(*(current_fp_nodal[lev][idim]),
                                                        current_fp_nodal[lev][idim]->nGrowVect(),
                                                        WarpX::do_single_precision_comms,
                                                        curr_period,
                                                        true);

        }
    }

    // set the boundary potentials appropriately
    setVectorPotentialBC(current_fp_nodal);

    // Compute the potential phi, by solving the Poisson equation
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( !IsPythonCallBackInstalled("poissonsolver"),
        "Python Level Poisson Solve not supported for Magnetostatic implementation.");

    computeVectorPotential( current_fp_nodal, vector_potential_fp_nodal, self_fields_required_precision,
                     self_fields_absolute_tolerance, self_fields_max_iters,
                     self_fields_verbosity );
}

/* Compute the vector potential `A` by solving the Poisson equation with `J` as
   a source.
   This uses the amrex solver.

   More specifically, this solves the equation
   \f[
       \vec{\nabla}^2 r \phi - (\vec{\beta}\cdot\vec{\nabla})^2 r \phi = -\frac{r \rho}{\epsilon_0}
   \f]

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
   \param[in] beta Represents the velocity of the source of `phi`
   \param[in] required_precision The relative convergence threshold for the MLMG solver
   \param[in] absolute_tolerance The absolute convergence threshold for the MLMG solver
   \param[in] max_iters The maximum number of iterations allowed for the MLMG solver
   \param[in] verbosity The verbosity setting for the MLMG solver
*/
void
WarpX::computeVectorPotential (const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>,3> >& curr,
                                amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>,3> >& A,
                                Real const required_precision,
                                Real absolute_tolerance,
                                int const max_iters,
                                int const verbosity) const
{
    // create a vector to our fields, sorted by level
    amrex::Vector<amrex::Array<amrex::MultiFab*,3>> sorted_curr;
    amrex::Vector<amrex::Array<amrex::MultiFab*,3>> sorted_A;
    for (int lev = 0; lev <= finest_level; ++lev) {
        sorted_curr.emplace_back(amrex::Array<amrex::MultiFab*,3> ({curr[lev][0].get(),
                                                                    curr[lev][1].get(),
                                                                    curr[lev][2].get()}));
        sorted_A.emplace_back(amrex::Array<amrex::MultiFab*,3> ({A[lev][0].get(),
                                                                 A[lev][1].get(),
                                                                 A[lev][2].get()}));
    }

#if defined(AMREX_USE_EB)
    std::optional<MagnetostaticSolver::EBCalcBfromVectorPotentialPerLevel> post_A_calculation({Efield_fp,
                                                                                               Bfield_fp,
                                                                                               vector_potential_grad_buf_e_stag,
                                                                                               vector_potential_grad_buf_b_stag,
                                                                                               vector_potential_fp_nodal,
                                                                                               vector_potential_old_fp_nodal,
                                                                                               dt[0]});

    amrex::Vector<amrex::EBFArrayBoxFactory const *> factories;
    for (int lev = 0; lev <= finest_level; ++lev) {
        factories.push_back(&WarpX::fieldEBFactory(lev));
    }
    std::optional<amrex::Vector<amrex::EBFArrayBoxFactory const *> > eb_farray_box_factory({factories});
#else
    std::optional<MagnetostaticSolver::EBCalcBfromVectorPotentialPerLevel> post_A_calculation;
    std::optional<amrex::Vector<amrex::FArrayBoxFactory const *> > eb_farray_box_factory;
#endif

    Magnetostatic::computeVectorPotential(
        sorted_curr,
        sorted_A,
        required_precision,
        absolute_tolerance,
        max_iters,
        verbosity,
        this->geom,
        this->dmap,
        this->grids,
        this->m_vector_poisson_boundary_handler,
        WarpX::do_single_precision_comms,
        this->ref_ratio,
        post_A_calculation,
        gett_new(0),
        eb_farray_box_factory
    );

}


/* \brief Set Dirichlet/Neumann boundary conditions for the magnetostatic solver.

    The given potential's values are fixed on the boundaries of the given
    dimension according to the desired values from the simulation input file,
    boundary.potential_lo and boundary.potential_hi.

   \param[inout] A The vector potential
*/
void
WarpX::setVectorPotentialBC ( amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>,3>>& A ) const
{
    // check if any dimension has non-periodic boundary conditions
    if (!m_vector_poisson_boundary_handler.has_non_periodic) return;

    auto dirichlet_flag = m_vector_poisson_boundary_handler.dirichlet_flag;

    // loop over all mesh refinement levels and set the boundary values
    for (int lev=0; lev <= max_level; lev++) {

        amrex::Box domain = Geom(lev).Domain();
        domain.surroundingNodes();
        for (int adim=0; adim < 3; adim++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(*A[lev][adim], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                // Extract the vector potential
                auto A_arr = A[lev][adim]->array(mfi);
                // Extract tileboxes for which to loop
                const Box& tb  = mfi.tilebox( A[lev][adim]->ixType().toIntVect() );

                // loop over dimensions
                for (int idim=0; idim<AMREX_SPACEDIM; idim++){
                    // check if neither boundaries in this dimension should be set
                    if (!(dirichlet_flag[adim][2*idim] || dirichlet_flag[adim][2*idim+1])) continue;

                    // a check can be added below to test if the boundary values
                    // are already correct, in which case the ParallelFor over the
                    // cells can be skipped

                    if (!domain.strictly_contains(tb)) {
                        amrex::ParallelFor( tb,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                                IntVect iv(AMREX_D_DECL(i,j,k));

                                if (dirichlet_flag[adim][2*idim] && iv[idim] == domain.smallEnd(idim))
                                    A_arr(i,j,k) = 0.;

                                if (dirichlet_flag[adim][2*idim+1] && iv[idim] == domain.bigEnd(idim))
                                    A_arr(i,j,k) = 0.;
                            } // loop ijk
                        );
                    }
                } // idim
            } // MFIter
        }
    } // lev
}

void MagnetostaticSolver::VectorPoissonBoundaryHandler::defineVectorPotentialBCs ( )
{
    for (int adim = 0; adim < 3; adim++) {
        int dim_start = 0;
#ifdef WARPX_DIM_RZ
        WarpX& warpx = WarpX::GetInstance();
        auto geom = warpx.Geom(0);
        if (geom.ProbLo(0) == 0){
            lobc[adim][0] = LinOpBCType::Neumann;
            dirichlet_flag[adim][0] = false;
            dim_start = 1;

            // handle the r_max boundary explicity
            if (WarpX::field_boundary_hi[0] == FieldBoundaryType::PEC) {
                if (adim == 0) {
                    hibc[adim][0] = LinOpBCType::Neumann;
                    dirichlet_flag[adim][1] = false;
                } else{
                    hibc[adim][0] = LinOpBCType::Dirichlet;
                    dirichlet_flag[adim][1] = true;
                }
            }
            else if (WarpX::field_boundary_hi[0] == FieldBoundaryType::Neumann) {
                hibc[adim][0] = LinOpBCType::Neumann;
                dirichlet_flag[adim][1] = false;
            }
        }
#endif
        for (int idim=dim_start; idim<AMREX_SPACEDIM; idim++){
            bool ndotA = (adim == idim);

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            if (idim == 1) ndotA = (adim == 2);
#endif

            if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::Periodic
                && WarpX::field_boundary_hi[idim] == FieldBoundaryType::Periodic ) {
                lobc[adim][idim] = LinOpBCType::Periodic;
                hibc[adim][idim] = LinOpBCType::Periodic;
                dirichlet_flag[adim][idim*2] = false;
                dirichlet_flag[adim][idim*2+1] = false;
            }
            else {
                has_non_periodic = true;
                if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC ) {
                    if (ndotA) {
                        lobc[adim][idim] = LinOpBCType::Neumann;
                        dirichlet_flag[adim][idim*2] = false;
                    } else {
                        lobc[adim][idim] = LinOpBCType::Dirichlet;
                        dirichlet_flag[adim][idim*2] = true;
                    }
                }

                else if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::Neumann ) {
                    lobc[adim][idim] = LinOpBCType::Neumann;
                    dirichlet_flag[adim][idim*2] = false;
                }
                else {
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                        "Field boundary conditions have to be either periodic, PEC, or neumann "
                        "when using the magnetostatic solver"
                    );
                }

                if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC ) {
                    if (ndotA) {
                        hibc[adim][idim] = LinOpBCType::Neumann;
                        dirichlet_flag[adim][idim*2+1] = false;
                    } else {
                        hibc[adim][idim] = LinOpBCType::Dirichlet;
                        dirichlet_flag[adim][idim*2+1] = true;
                    }
                }
                else if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::Neumann ) {
                    hibc[adim][idim] = LinOpBCType::Neumann;
                    dirichlet_flag[adim][idim*2+1] = false;
                }
                else {
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                        "Field boundary conditions have to be either periodic, PEC, or neumann "
                        "when using the electrostatic solver"
                    );
                }
            }
        }
    }
    bcs_set = true;
}


void MagnetostaticSolver::EBCalcBfromVectorPotentialPerLevel::doInterp(const int isrc, const int idst, const int lev)
{

    WarpX &warpx = WarpX::GetInstance();

    // Grab Interpolation Coefficients
    // Order of finite-order centering of fields
    const int fg_nox = WarpX::field_centering_nox;
    const int fg_noy = WarpX::field_centering_noy;
    const int fg_noz = WarpX::field_centering_noz;

    // Device vectors of stencil coefficients used for finite-order centering of fields
    amrex::Real const * stencil_coeffs_x = warpx.device_field_centering_stencil_coeffs_x.data();
    amrex::Real const * stencil_coeffs_y = warpx.device_field_centering_stencil_coeffs_y.data();
    amrex::Real const * stencil_coeffs_z = warpx.device_field_centering_stencil_coeffs_z.data();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*(m_grad_buf_b_stag[lev][idst]), TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        IntVect const src_stag = m_grad_buf_e_stag[lev][isrc]->ixType().toIntVect();
        IntVect const dst_stag = m_grad_buf_b_stag[lev][idst]->ixType().toIntVect();

        Array4<amrex::Real const> const& src_arr = m_grad_buf_e_stag[lev][isrc]->const_array(mfi);
        Array4<amrex::Real> const& dst_arr = m_grad_buf_b_stag[lev][idst]->array(mfi);

        // Loop includes ghost cells (`growntilebox`)
        // (input arrays will be padded with zeros beyond ghost cells
        // for out-of-bound accesses due to large-stencil operations)
        Box bx = mfi.tilebox(dst_stag, m_grad_buf_b_stag[lev][idst]->nGrowVect());

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
        {
            warpx_interp(j, k, l, dst_arr, src_arr, dst_stag, src_stag, fg_nox, fg_noy, fg_noz,
                stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
        });
    }
}

void MagnetostaticSolver::EBCalcBfromVectorPotentialPerLevel::doEfieldCalc(const int lev)
{
    amrex::Real c = PhysConst::c;
    amrex::Real cdti = 1_rt/(c*m_dt);

    for (int adim=0; adim<3; adim++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*(m_e_field[lev][adim]), TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<amrex::Real> const& Efield_arr = m_e_field[lev][adim]->array(mfi);
            Array4<amrex::Real const> const& A_arr = m_A[lev][adim]->const_array(mfi);
            Array4<amrex::Real> const& A_old_arr = m_A_old[lev][adim]->array(mfi);

            IntVect const stag = m_e_field[lev][adim]->ixType().toIntVect();

            Box bx = mfi.tilebox(stag);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
            {
                // Calculate E-field contribution (-1/c dA/dt)
                Efield_arr(j,k,l) -= (A_arr(j,k,l)-A_old_arr(j,k,l))*cdti;

                // Update old buffer
                A_old_arr(j,k,l) = A_arr(j,k,l);
            });
        }
    }
}

void MagnetostaticSolver::EBCalcBfromVectorPotentialPerLevel::operator()(amrex::Array<std::unique_ptr<amrex::MLMG>,3> & mlmg, int const lev)
{
    using namespace amrex::literals;

    // This operator gets the gradient solution on the cell edges, aligned with E field staggered grid
    // This routine interpolates to the B-field staggered grid,

    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> buf_ptr =
    {
#if defined(WARPX_DIM_3D)
        m_grad_buf_e_stag[lev][0].get(),
        m_grad_buf_e_stag[lev][1].get(),
        m_grad_buf_e_stag[lev][2].get()
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        m_grad_buf_e_stag[lev][0].get(),
        m_grad_buf_e_stag[lev][2].get()
#endif
    };

    // This will grab the gradient values for Ax
    mlmg[0]->getGradSolution({buf_ptr});

    // Interpolate dAx/dz to By grid buffer, then add to By
    this->doInterp(2,1,lev);
    MultiFab::Add(*(m_b_field[lev][1]), *(m_grad_buf_b_stag[lev][1]), 0, 0, 1, m_b_field[lev][1]->nGrowVect() );

    // Interpolate dAx/dy to Bz grid buffer, then subtract from Bz
    this->doInterp(1,2,lev);
    m_grad_buf_b_stag[lev][2]->mult(-1._rt);
    MultiFab::Add(*(m_b_field[lev][2]), *(m_grad_buf_b_stag[lev][2]), 0, 0, 1, m_b_field[lev][2]->nGrowVect() );

    // This will grab the gradient values for Ay
    mlmg[1]->getGradSolution({buf_ptr});

    // Interpolate dAy/dx to Bz grid buffer, then add to Bz
    this->doInterp(0,2,lev);
    MultiFab::Add(*(m_b_field[lev][2]), *(m_grad_buf_b_stag[lev][2]), 0, 0, 1, m_b_field[lev][2]->nGrowVect() );

    // Interpolate dAy/dz to Bx grid buffer, then subtract from Bx
    this->doInterp(2,0,lev);
    m_grad_buf_b_stag[lev][0]->mult(-1._rt);
    MultiFab::Add(*(m_b_field[lev][0]), *(m_grad_buf_b_stag[lev][0]), 0, 0, 1, m_b_field[lev][0]->nGrowVect() );

    // This will grab the gradient values for Az
    mlmg[2]->getGradSolution({buf_ptr});

    // Interpolate dAz/dy to Bx grid buffer, then add to Bx
    this->doInterp(1,0,lev);
    MultiFab::Add(*(m_b_field[lev][0]), *(m_grad_buf_b_stag[lev][0]), 0, 0, 1, m_b_field[lev][0]->nGrowVect() );

    // Interpolate dAz/dx to By grid buffer, then subtract from By
    this->doInterp(0,1,lev);
    m_grad_buf_b_stag[lev][1]->mult(-1._rt);
    MultiFab::Add(*(m_b_field[lev][1]), *(m_grad_buf_b_stag[lev][1]), 0, 0, 1, m_b_field[lev][1]->nGrowVect() );

    // Additionally compute the -1/c dA/dt term and add to Efield
    this->doEfieldCalc(lev);
}