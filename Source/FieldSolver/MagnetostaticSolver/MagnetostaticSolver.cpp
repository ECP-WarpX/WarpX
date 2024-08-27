/* Copyright 2022-2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: S. Eric Clark (LLNL), Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "FieldSolver/MagnetostaticSolver/MagnetostaticSolver.H"
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Python/callbacks.H"
#include "Utils/WarpXConst.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Parallelization/WarpXComm_K.H"

#include <ablastr/utils/Communication.H>
#include <ablastr/warn_manager/WarnManager.H>
#include <ablastr/fields/VectorPoissonSolver.H>

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
WarpX::ComputeMagnetostaticField()
{
    WARPX_PROFILE("WarpX::ComputeMagnetostaticField");
    // Fields have been reset in Electrostatic solver for this time step, these fields
    // are added into the B fields after electrostatic solve

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        this->max_level == 0,
        "Magnetostatic solver not implemented with mesh refinement."
    );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        do_magnetostatic_solve,
        "Magnetostatic solver called but warpx.do_magnetostatic is false."
    );

    ComputeMagnetostaticFieldLabFrame();
}

void
WarpX::ComputeMagnetostaticFieldLabFrame()
{
    WARPX_PROFILE("WarpX::ComputeMagnetostaticFieldLabFrame");

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(n_rz_azimuthal_modes == 1,
                                     "Error: RZ magnetostatic only implemented for a single mode");
#endif

    // Perform current deposition at t_{n+1/2} (source of Poisson solver).
    mypc->DepositCurrent(current_fp, dt[0], -0.5_rt * dt[0]);

    // Synchronize J and rho:
    // filter (if used), exchange guard cells, interpolate across MR levels
    // and apply boundary conditions
    SyncCurrent(current_fp, current_cp, current_buf);
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        ApplyJfieldBoundary(
            lev, current_fp[lev][0].get(), current_fp[lev][1].get(),
            current_fp[lev][2].get(), PatchType::fine
        );
        if (lev > 0) {
            ApplyJfieldBoundary(
                lev, current_cp[lev][0].get(), current_cp[lev][1].get(),
                current_cp[lev][2].get(), PatchType::coarse
            );
        }
    }

    // SyncCurrent does not include a call to FillBoundary
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            current_fp[lev][idim]->FillBoundary(Geom(lev).periodicity());
        }
    }

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!m_vector_poisson_boundary_handler.bcs_set) {
        m_vector_poisson_boundary_handler.defineVectorPotentialBCs();
    }
    // set the boundary and current density potentials
    setVectorPotentialBC(Afield_fp_nodal);

    // Compute the vector potential A, by solving the Poisson equation
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( !IsPythonCallbackInstalled("poissonsolver"),
        "Python Level Poisson Solve not supported for Magnetostatic implementation.");

    const amrex::Real magnetostatic_absolute_tolerance = self_fields_absolute_tolerance*PhysConst::c;

    computeVectorPotential(
        current_fp, Afield_fp_nodal,
        self_fields_required_precision, magnetostatic_absolute_tolerance,
        self_fields_max_iters, self_fields_verbosity
    );

    // Grab Interpolation Coefficients
    // Order of finite-order centering of fields
    const int fg_nox = WarpX::field_centering_nox;
    const int fg_noy = WarpX::field_centering_noy;
    const int fg_noz = WarpX::field_centering_noz;

    // Device vectors of stencil coefficients used for finite-order centering of fields
    amrex::Real const * stencil_coeffs_x = WarpX::device_field_centering_stencil_coeffs_x.data();
    amrex::Real const * stencil_coeffs_y = WarpX::device_field_centering_stencil_coeffs_y.data();
    amrex::Real const * stencil_coeffs_z = WarpX::device_field_centering_stencil_coeffs_z.data();

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // interpolate Afield_fp_nodal[lev][i] to Afield_fp[lev][i]
        for (int i = 0; i < 3; ++i) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*Afield_fp[lev][i], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                IntVect const src_stag = Afield_fp_nodal[lev][i]->ixType().toIntVect();
                IntVect const dst_stag = Afield_fp[lev][i]->ixType().toIntVect();

                Array4<amrex::Real const> const& src_arr = Afield_fp_nodal[lev][i]->const_array(mfi);
                Array4<amrex::Real> const& dst_arr = Afield_fp[lev][i]->array(mfi);

                const Box bx = mfi.tilebox();

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, dst_arr, src_arr, dst_stag, src_stag, fg_nox, fg_noy, fg_noz,
                        stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
                });
            }
        }
    }

    ComputeBfromVectorPotential();
}

/* Compute the vector potential `A` by solving the Poisson equation with `J` as
   a source.
   This uses the amrex solver.

    More specifically, this solves the equation
    \f[
        \vec{\nabla}^2 r \vec{A} = - r \mu_0 \vec{J}
 \f]

   \param[in] curr The current density
   \param[out] A The vector potential to be computed by this function
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
    amrex::Vector<amrex::EBFArrayBoxFactory const *> factories;
    for (int lev = 0; lev <= finest_level; ++lev) {
        factories.push_back(&WarpX::fieldEBFactory(lev));
    }
    const std::optional<amrex::Vector<amrex::EBFArrayBoxFactory const *> > eb_farray_box_factory({factories});
#else
    const std::optional<amrex::Vector<amrex::FArrayBoxFactory const *> > eb_farray_box_factory;
#endif

    ablastr::fields::computeVectorPotential(
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
        std::nullopt, // post_A_calculation,
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
    if (!m_vector_poisson_boundary_handler.has_non_periodic) { return; }

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
                const Box& tb  = mfi.tilebox( A[lev][adim]->ixType().toIntVect());

                // loop over dimensions
                for (int idim=0; idim<AMREX_SPACEDIM; idim++){
                    // check if neither boundaries in this dimension should be set
                    if (!(dirichlet_flag[adim][2*idim] || dirichlet_flag[adim][2*idim+1])) { continue; }

                    // a check can be added below to test if the boundary values
                    // are already correct, in which case the ParallelFor over the
                    // cells can be skipped

                    if (!domain.strictly_contains(tb)) {
                        amrex::ParallelFor( tb,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                                IntVect iv(AMREX_D_DECL(i,j,k));

                                if (dirichlet_flag[adim][2*idim] && iv[idim] == domain.smallEnd(idim)) {
                                    A_arr(i,j,k) = 0.;
                                }

                                if (dirichlet_flag[adim][2*idim+1] && iv[idim] == domain.bigEnd(idim)) {
                                    A_arr(i,j,k) = 0.;
                                }
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
#ifdef WARPX_DIM_RZ
        int dim_start = 0;
        WarpX& warpx = WarpX::GetInstance();
        auto geom = warpx.Geom(0);
        if (geom.ProbLo(0) == 0){
            lobc[adim][0] = LinOpBCType::Neumann;
            dirichlet_flag[adim][0] = false;
            dim_start = 1;

            // handle the r_max boundary explicitly
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
#else
        const int dim_start = 0;
#endif
        bool ndotA = false;
        for (int idim=dim_start; idim<AMREX_SPACEDIM; idim++){
            ndotA = (adim == idim);

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            if (idim == 1) { ndotA = (adim == 2); }
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
                    WARPX_ABORT_WITH_MESSAGE(
                        "Field boundary conditions have to be either periodic, PEC or neumann "
                        "when using the magnetostatic solver,  but they are " + GetFieldBCTypeString(WarpX::field_boundary_lo[idim])
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
                    WARPX_ABORT_WITH_MESSAGE(
                        "Field boundary conditions have to be either periodic, PEC or neumann "
                        "when using the magnetostatic solver,  but they are " + GetFieldBCTypeString(WarpX::field_boundary_lo[idim])
                    );
                }
            }
        }
    }
    bcs_set = true;
}
