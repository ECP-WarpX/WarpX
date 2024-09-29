/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
/*
 * This file was copied and edited from PoissonSolver.H in the same directory.
 */
#ifndef ABLASTR_SEMI_IMPLICIT_POISSON_SOLVER_H
#define ABLASTR_SEMI_IMPLICIT_POISSON_SOLVER_H

#include <ablastr/constant.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/utils/TextMsg.H>
#include <ablastr/warn_manager/WarnManager.H>
#include <ablastr/math/fft/AnyFFT.H>
#include <ablastr/fields/Interpolate.H>
#include <ablastr/profiler/ProfilerWrapper.H>

#if defined(WARPX_USE_FFT) && defined(WARPX_DIM_3D)
#include <ablastr/fields/IntegratedGreenFunctionSolver.H>
#endif

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
#include <AMReX_MFInterp_C.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>
#include <AMReX_MLEBNodeFDLaplacian.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EBFabFactory.H>
#endif

#include <array>
#include <optional>


namespace ablastr::fields {

/** Compute the potential `phi` by solving the Poisson equation with a modifed dielectric function
 *
 * Uses `rho` as a source. This uses the AMReX solver.
 *
 * More specifically, this solves the equation
 * \f[
 *   \nabla \cdot \sigma \nabla \phi = - \rho/\epsilon_0
 * \f]
 *
 * \tparam T_BoundaryHandler handler for boundary conditions, for example @see ElectrostaticSolver::PoissonBoundaryHandler
 * \tparam T_PostPhiCalculationFunctor a calculation per level directly after phi was calculated
 * \tparam T_FArrayBoxFactory usually nothing or an amrex::EBFArrayBoxFactory (EB ONLY)
 * \param[in] rho The charge density a given species
 * \param[out] phi The potential to be computed by this function
 * \param[in] sigma The semi-implicit matrix
 * \param[in] relative_tolerance The relative convergence threshold for the MLMG solver
 * \param[in] absolute_tolerance The absolute convergence threshold for the MLMG solver
 * \param[in] max_iters The maximum number of iterations allowed for the MLMG solver
 * \param[in] verbosity The verbosity setting for the MLMG solver
 * \param[in] geom the geometry per level (e.g., from AmrMesh)
 * \param[in] dmap the distribution mapping per level (e.g., from AmrMesh)
 * \param[in] grids the grids per level (e.g., from AmrMesh)
 * \param[in] boundary_handler a handler for boundary conditions, for example @see ElectrostaticSolver::PoissonBoundaryHandler
 * \param[in] is_solver_igf_on_lev0 boolean to select the Poisson solver: 1 for FFT on level 0 & Multigrid on other levels, 0 for Multigrid on all levels
 * \param[in] do_single_precision_comms perform communications in single precision
 * \param[in] rel_ref_ratio mesh refinement ratio between levels (default: 1)
 * \param[in] post_phi_calculation perform a calculation per level directly after phi was calculated; required for embedded boundaries (default: none)
 * \param[in] current_time the current time; required for embedded boundaries (default: none)
 * \param[in] eb_farray_box_factory a factory for field data, @see amrex::EBFArrayBoxFactory; required for embedded boundaries (default: none)
 */
template<
    typename T_BoundaryHandler,
    typename T_PostPhiCalculationFunctor = std::nullopt_t,
    typename T_FArrayBoxFactory = void
>
void
computeSemiImplicitPhi (
    ablastr::fields::MultiLevelScalarField const& rho,
    ablastr::fields::MultiLevelScalarField const& phi,
    amrex::MultiFab const & sigma,
    amrex::Real relative_tolerance,
    amrex::Real absolute_tolerance,
    int max_iters,
    int verbosity,
    amrex::Vector<amrex::Geometry> const& geom,
    amrex::Vector<amrex::DistributionMapping> const& dmap,
    amrex::Vector<amrex::BoxArray> const& grids,
    [[maybe_unused]] utils::enums::GridType grid_type,
    T_BoundaryHandler const boundary_handler,
    bool is_solver_igf_on_lev0,
    [[maybe_unused]] bool eb_enabled = false,
    bool do_single_precision_comms = false,
    std::optional<amrex::Vector<amrex::IntVect> > rel_ref_ratio = std::nullopt,
    [[maybe_unused]] T_PostPhiCalculationFunctor post_phi_calculation = std::nullopt,
    [[maybe_unused]] std::optional<amrex::Real const> current_time = std::nullopt, // only used for EB
    [[maybe_unused]] std::optional<amrex::Vector<T_FArrayBoxFactory const *> > eb_farray_box_factory = std::nullopt // only used for EB
) {
    using namespace amrex::literals;

    ABLASTR_PROFILE("computeSemiImplicitPhi");

    if (!rel_ref_ratio.has_value()) {
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(rho.size() == 1u,
                                           "rel_ref_ratio must be set if mesh-refinement is used");
        rel_ref_ratio = amrex::Vector<amrex::IntVect>{{amrex::IntVect(AMREX_D_DECL(1, 1, 1))}};
    }
    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE( !is_solver_igf_on_lev0,
        "FFT solver cannot be used with semi-implicit Poisson solve");

    auto const finest_level = static_cast<int>(rho.size() - 1);

    // determine if rho is zero everywhere
    amrex::Real max_norm_b = 0.0;
    for (int lev=0; lev<=finest_level; lev++) {
        max_norm_b = amrex::max(max_norm_b, rho[lev]->norm0());
    }
    amrex::ParallelDescriptor::ReduceRealMax(max_norm_b);

    const bool always_use_bnorm = (max_norm_b > 0);
    if (!always_use_bnorm) {
        if (absolute_tolerance == 0.0) { absolute_tolerance = amrex::Real(1e-6); }
        ablastr::warn_manager::WMRecordWarning(
                "ElectrostaticSolver",
                "Max norm of rho is 0",
                ablastr::warn_manager::WarnPriority::low
        );
    }

    const amrex::LPInfo info;

    for (int lev=0; lev<=finest_level; lev++) {

        // Use the Multigrid (MLMG) solver but first scale rho appropriately
        using namespace ablastr::constant::SI;
        rho[lev]->mult(-1._rt/ep0);  // TODO: when do we "un-multiply" this? We need to document this side-effect!

#if defined(AMREX_USE_EB) || defined(WARPX_DIM_RZ)
        // Use the EB enabled solver.
        amrex::MLEBNodeFDLaplacian linop( {geom[lev]}, {grids[lev]}, {dmap[lev]}, info
#if defined(AMREX_USE_EB)
            , {eb_farray_box_factory.value()[lev]}
#endif
        );
#else
        amrex::MLNodeLaplacian linop( {geom[lev]}, {grids[lev]}, {dmap[lev]}, info);
#endif

#if defined(AMREX_USE_EB)
        // if the EB potential only depends on time, the potential can be passed
        // as a float instead of a callable
        if (boundary_handler.phi_EB_only_t) {
            linop.setEBDirichlet(boundary_handler.potential_eb_t(current_time.value()));
        }
        else
            linop.setEBDirichlet(boundary_handler.getPhiEB(current_time.value()));
#endif
        linop.setSigma(lev, sigma);

        // Solve the Poisson equation
        linop.setDomainBC( boundary_handler.lobc, boundary_handler.hibc );
#ifdef WARPX_DIM_RZ
        linop.setRZ(true);
#endif
        amrex::MLMG mlmg(linop); // actual solver defined here
        mlmg.setVerbose(verbosity);
        mlmg.setMaxIter(max_iters);
        mlmg.setAlwaysUseBNorm(always_use_bnorm);

        // Solve Poisson equation at lev
        mlmg.solve( {phi[lev]}, {rho[lev]},
                    relative_tolerance, absolute_tolerance );

        // needed for solving the levels by levels:
        // - coarser level is initial guess for finer level
        // - coarser level provides boundary values for finer level patch
        // Interpolation from phi[lev] to phi[lev+1]
        // (This provides both the boundary conditions and initial guess for phi[lev+1])
        if (lev < finest_level) {

            // Allocate phi_cp for lev+1
            amrex::BoxArray ba = phi[lev+1]->boxArray();
            const amrex::IntVect& refratio = rel_ref_ratio.value()[lev];
            ba.coarsen(refratio);
            const int ncomp = linop.getNComp();
            amrex::MultiFab phi_cp(ba, phi[lev+1]->DistributionMap(), ncomp, 1);

            // Copy from phi[lev] to phi_cp (in parallel)
            const amrex::IntVect& ng = amrex::IntVect::TheUnitVector();
            const amrex::Periodicity& crse_period = geom[lev].periodicity();

            ablastr::utils::communication::ParallelCopy(
                phi_cp,
                *phi[lev],
                0,
                0,
                1,
                ng,
                ng,
                do_single_precision_comms,
                crse_period
            );

            // Local interpolation from phi_cp to phi[lev+1]
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(*phi[lev + 1], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Array4<amrex::Real> const phi_fp_arr = phi[lev + 1]->array(mfi);
                amrex::Array4<amrex::Real const> const phi_cp_arr = phi_cp.array(mfi);

                details::PoissonInterpCPtoFP const interp(phi_fp_arr, phi_cp_arr, refratio);

                amrex::Box const b = mfi.tilebox(phi[lev + 1]->ixType().toIntVect());
                amrex::ParallelFor(b, interp);
            }

        }

        // Run additional operations, such as calculation of the E field for embedded boundaries
        if constexpr (!std::is_same<T_PostPhiCalculationFunctor, std::nullopt_t>::value) {
            if (post_phi_calculation.has_value()) {
                post_phi_calculation.value()(mlmg, lev);
            }
        }

    } // loop over lev(els)
}

} // namespace ablastr::fields

#endif // ABLASTR_SEMI_IMPLICIT_POISSON_SOLVER_H
