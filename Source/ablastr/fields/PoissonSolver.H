/* Copyright 2019-2022 Axel Huebl, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_POISSON_SOLVER_H
#define ABLASTR_POISSON_SOLVER_H

#include "Utils/WarpXConst.H"

#include <ablastr/utils/Communication.H>
#include <ablastr/utils/TextMsg.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_MLNodeTensorLaplacian.H>
#if defined(AMREX_USE_EB) || defined(WARPX_DIM_RZ)
#   include <AMReX_MLEBNodeFDLaplacian.H>
#endif
#ifdef AMREX_USE_EB
#   include <AMReX_EBFabFactory.H>
#endif
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_MFInterp_C.H>

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

#include <array>
#include <optional>


namespace ablastr::fields {

namespace details
{
    /** Local interpolation from phi_cp to phi[lev+1]
     *
     * This is needed to work-around an NVCC limitation in downstream code (ImpactX),
     * when nesting lambdas. Otherwise this could be written directly into the
     * ParallelFor.
     *
     * @param[out] phi_fp_arr phi on the fine level
     * @param[in] phi_cp_arr phi on the coarse level
     * @param[in] refratio refinement ration
     */
    struct PoissonInterpCPtoFP
    {
        PoissonInterpCPtoFP(
                amrex::Array4<amrex::Real> const phi_fp_arr,
                amrex::Array4<amrex::Real const> const phi_cp_arr,
                amrex::IntVect const refratio)
        : m_phi_fp_arr(phi_fp_arr), m_phi_cp_arr(phi_cp_arr), m_refratio(refratio)
        {}

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void
        operator() (long i, long j, long k) const noexcept
        {
            amrex::mf_nodebilin_interp(i, j, k, 0, m_phi_fp_arr, 0, m_phi_cp_arr,
                                       0, m_refratio);
        }

        amrex::Array4<amrex::Real> const m_phi_fp_arr;
        amrex::Array4<amrex::Real const> const m_phi_cp_arr;
        amrex::IntVect const m_refratio;
    };
}

/** Compute the potential `phi` by solving the Poisson equation
 *
 * Uses `rho` as a source, assuming that the source moves at a
 * constant speed \f$\vec{\beta}\f$. This uses the AMReX solver.
 *
 * More specifically, this solves the equation
 * \f[
 *   \vec{\nabla}^2 r \phi - (\vec{\beta}\cdot\vec{\nabla})^2 r \phi = -\frac{r \rho}{\epsilon_0}
 * \f]
 *
 * \tparam T_BoundaryHandler handler for boundary conditions, for example @see ElectrostaticSolver::PoissonBoundaryHandler
 * \tparam T_PostPhiCalculationFunctor a calculation per level directly after phi was calculated
 * \tparam T_FArrayBoxFactory usually nothing or an amrex::EBFArrayBoxFactory (EB ONLY)
 * \param[in] rho The charge density a given species
 * \param[out] phi The potential to be computed by this function
 * \param[in] beta Represents the velocity of the source of `phi`
 * \param[in] relative_tolerance The relative convergence threshold for the MLMG solver
 * \param[in] absolute_tolerance The absolute convergence threshold for the MLMG solver
 * \param[in] max_iters The maximum number of iterations allowed for the MLMG solver
 * \param[in] verbosity The verbosity setting for the MLMG solver
 * \param[in] geom the geometry per level (e.g., from AmrMesh)
 * \param[in] dmap the distribution mapping per level (e.g., from AmrMesh)
 * \param[in] grids the grids per level (e.g., from AmrMesh)
 * \param[in] boundary_handler a handler for boundary conditions, for example @see ElectrostaticSolver::PoissonBoundaryHandler
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
computePhi (amrex::Vector<amrex::MultiFab*> const & rho,
            amrex::Vector<amrex::MultiFab*> & phi,
            std::array<amrex::Real, 3> const beta,
            amrex::Real const relative_tolerance,
            amrex::Real absolute_tolerance,
            int const max_iters,
            int const verbosity,
            amrex::Vector<amrex::Geometry> const geom,
            amrex::Vector<amrex::DistributionMapping> const dmap,
            amrex::Vector<amrex::BoxArray> const grids,
            T_BoundaryHandler const boundary_handler,
            bool const do_single_precision_comms = false,
            std::optional<amrex::Vector<amrex::IntVect> > rel_ref_ratio = std::nullopt,
            [[maybe_unused]] T_PostPhiCalculationFunctor post_phi_calculation = std::nullopt,
            [[maybe_unused]] std::optional<amrex::Real const> current_time = std::nullopt, // only used for EB
            [[maybe_unused]] std::optional<amrex::Vector<T_FArrayBoxFactory const *> > eb_farray_box_factory = std::nullopt // only used for EB
)
{
    using namespace amrex::literals;

    if (!rel_ref_ratio.has_value()) {
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(rho.size() == 1u,
                                           "rel_ref_ratio must be set if mesh-refinement is used");
        rel_ref_ratio = amrex::Vector<amrex::IntVect>{{amrex::IntVect(AMREX_D_DECL(1, 1, 1))}};
    }

    int const finest_level = rho.size() - 1u;

    // scale rho appropriately; also determine if rho is zero everywhere
    amrex::Real max_norm_b = 0.0;
    for (int lev=0; lev<=finest_level; lev++) {
        rho[lev]->mult(-1._rt/PhysConst::ep0); // TODO: when do we "un-multiply" this? We need to document this side-effect!
        max_norm_b = amrex::max(max_norm_b, rho[lev]->norm0());
    }
    amrex::ParallelDescriptor::ReduceRealMax(max_norm_b);

    bool always_use_bnorm = (max_norm_b > 0);
    if (!always_use_bnorm) {
        if (absolute_tolerance == 0.0) absolute_tolerance = amrex::Real(1e-6);
        ablastr::warn_manager::WMRecordWarning(
                "ElectrostaticSolver",
                "Max norm of rho is 0",
                ablastr::warn_manager::WarnPriority::low
        );
    }

    amrex::LPInfo info;
    for (int lev=0; lev<=finest_level; lev++) {
        // Set the value of beta
        amrex::Array<amrex::Real,AMREX_SPACEDIM> beta_solver =
#if defined(WARPX_DIM_1D_Z)
                {{ beta[2] }};  // beta_x and beta_z
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                {{ beta[0], beta[2] }};  // beta_x and beta_z
#else
                {{ beta[0], beta[1], beta[2] }};
#endif

#if !(defined(AMREX_USE_EB) || defined(WARPX_DIM_RZ))
        // Determine whether to use semi-coarsening
        amrex::Array<amrex::Real,AMREX_SPACEDIM> dx_scaled
                {AMREX_D_DECL(geom[lev].CellSize(0)/std::sqrt(1._rt-beta_solver[0]*beta_solver[0]),
                              geom[lev].CellSize(1)/std::sqrt(1._rt-beta_solver[1]*beta_solver[1]),
                              geom[lev].CellSize(2)/std::sqrt(1._rt-beta_solver[2]*beta_solver[2]))};
        int max_semicoarsening_level = 0;
        int semicoarsening_direction = -1;
        int min_dir = std::distance(dx_scaled.begin(),
                                    std::min_element(dx_scaled.begin(),dx_scaled.end()));
        int max_dir = std::distance(dx_scaled.begin(),
                                    std::max_element(dx_scaled.begin(),dx_scaled.end()));
        if (dx_scaled[max_dir] > dx_scaled[min_dir]) {
            semicoarsening_direction = max_dir;
            max_semicoarsening_level = static_cast<int>
            (std::log2(dx_scaled[max_dir]/dx_scaled[min_dir]));
        }
        if (max_semicoarsening_level > 0) {
            info.setSemicoarsening(true);
            info.setMaxSemicoarseningLevel(max_semicoarsening_level);
            info.setSemicoarseningDirection(semicoarsening_direction);
        }
#endif

#if defined(AMREX_USE_EB) || defined(WARPX_DIM_RZ)
        // In the presence of EB or RZ: the solver assumes that the beam is
        // propagating along  one of the axes of the grid, i.e. that only *one*
        // of the components of `beta` is non-negligible.
        amrex::MLEBNodeFDLaplacian linop( {geom[lev]}, {grids[lev]}, {dmap[lev]}, info
#if defined(AMREX_USE_EB)
            , {eb_farray_box_factory.value()[lev]}
#endif
        );

        // Note: this assumes that the beam is propagating along
        // one of the axes of the grid, i.e. that only *one* of the
        // components of `beta` is non-negligible. // we use this
#if defined(WARPX_DIM_RZ)
        linop.setSigma({0._rt, 1._rt-beta_solver[1]*beta_solver[1]});
#else
        linop.setSigma({AMREX_D_DECL(
            1._rt-beta_solver[0]*beta_solver[0],
            1._rt-beta_solver[1]*beta_solver[1],
            1._rt-beta_solver[2]*beta_solver[2])});
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
#else
        // In the absence of EB and RZ: use a more generic solver
        // that can handle beams propagating in any direction
        amrex::MLNodeTensorLaplacian linop( {geom[lev]}, {grids[lev]},
                                     {dmap[lev]}, info );
        linop.setBeta( beta_solver ); // for the non-axis-aligned solver
#endif

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
        if constexpr (!std::is_same<T_PostPhiCalculationFunctor, std::nullopt_t>::value)
            if (post_phi_calculation.has_value())
                post_phi_calculation.value()(mlmg, lev);

    } // loop over lev(els)
}

} // namespace ablastr::fields

#endif // ABLASTR_POISSON_SOLVER_H
