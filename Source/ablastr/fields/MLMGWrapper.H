/* Copyright 2026 The WarpX Community
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_MLMG_WRAPPER_H_
#define ABLASTR_MLMG_WRAPPER_H_

#include <ablastr/fields/MLMGOptions.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_IntVect.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EBFabFactory.H>
#endif

#include <algorithm>
#include <array>
#include <cfenv>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#if defined(AMREX_USE_EB)

namespace ablastr::fields
{

/** Name of a bottom solver, for messages */
inline std::string bottomSolverName (amrex::BottomSolver s)
{
    switch (s) {
    case amrex::BottomSolver::smoother: return "smoother";
    case amrex::BottomSolver::bicgstab: return "bicgstab";
    case amrex::BottomSolver::cg:       return "cg";
    case amrex::BottomSolver::bicgcg:   return "bicgcg";
    case amrex::BottomSolver::cgbicg:   return "cgbicg";
    case amrex::BottomSolver::hypre:    return "hypre";
    case amrex::BottomSolver::petsc:    return "petsc";
    case amrex::BottomSolver::custom:   return "custom";
    case amrex::BottomSolver::algmg:    return "algmg";
    default:                            return "default";
    }
}

/** Operator and solver of one level's Poisson solve with an embedded boundary
 *
 * Holds an amrex::MLEBNodeFDLaplacian and its amrex::MLMG for one level in
 * Cartesian geometry.  Kept between solves by the caller, so that MLMG keeps
 * its setup (the AlgMG hierarchy in particular).
 *
 * With MLMGOptions::tune_mlmg, the constructor times trial solves over the
 * maximum coarsening level, the agglomeration grid size, the bottom solver
 * (the configured one, default, algmg, smoother and, on CUDA and HIP, custom), a
 * single-box layout and, on GPUs, the stream synchronization mode, and keeps
 * the fastest operator and solver.  Every configuration gets one untimed
 * warm-up solve, so that its setup cost does not count: the solver is reused
 * and the setup is paid once.
 */
class MLMGWrapper
{
public:
    /** The MLMG parameters that the tuning chooses from */
    struct Parms {
        bool single_box = false; //!< solve on one box covering the whole domain
        int mcl = 30;            //!< maximum coarsening level (AMReX default)
        int ags = -1;            //!< agglomeration grid size, -1 for the AMReX default
        bool no_gpu_sync = true; //!< GPU only: single stream, no stream synchronization
        amrex::BottomSolver bottom_solver = amrex::BottomSolver::Default;
    };

    using BCArray = amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>;
    using RealArray = amrex::Array<amrex::Real, AMREX_SPACEDIM>;

    /** Build, and tune if enabled, the operator and solver of one level
     *
     * \param[in] lev mesh refinement level
     * \param[in] geom, grids, dmap, factory the level's geometry, cell-centered grids, EB factory
     * \param[in] info the coarsening options; tuning may change the coarsening depth
     * \param[in] lobc, hibc the domain boundary conditions
     * \param[in] sigma the operator coefficients, (1 - beta^2) per direction
     * \param[in] options the MLMG options
     * \param[in] ng 1 if phi has ghost nodes that MLMG must fill
     * \param[in] set_eb_dirichlet callable setting the EB potential on an operator
     * \param[in] rho right-hand side of the trial solves (already scaled); nullptr: no tuning
     * \param[in] relative_tolerance, absolute_tolerance the tolerances of the trial solves
     */
    template <typename F>
    MLMGWrapper (int lev,
                 amrex::Geometry const& geom,
                 amrex::BoxArray const& grids,
                 amrex::DistributionMapping const& dmap,
                 amrex::EBFArrayBoxFactory const& factory,
                 amrex::LPInfo const& info,
                 BCArray const& lobc,
                 BCArray const& hibc,
                 RealArray const& sigma,
                 MLMGOptions const& options,
                 int ng,
                 F const& set_eb_dirichlet,
                 amrex::MultiFab const* rho,
                 amrex::Real relative_tolerance,
                 amrex::Real absolute_tolerance)
        : m_geom(geom), m_ng(ng)
    {
        m_parms.mcl = info.max_coarsening_level;
        m_parms.ags = info.agg_grid_size;
        m_parms.bottom_solver = options.bottom_solver;

        // Operator and solver on a layout, with given parameters
        auto const make = [&] (amrex::BoxArray const& ba, amrex::DistributionMapping const& dm,
                               amrex::EBFArrayBoxFactory const* fac, Parms const& parms)
        {
            amrex::LPInfo linfo = info;
            linfo.setMaxCoarseningLevel(parms.mcl).setAgglomerationGridSize(parms.ags);
            auto linop = std::make_unique<amrex::MLEBNodeFDLaplacian>();
            linop->define(amrex::Vector<amrex::Geometry>{geom}, amrex::Vector<amrex::BoxArray>{ba},
                          amrex::Vector<amrex::DistributionMapping>{dm}, linfo,
                          amrex::Vector<amrex::EBFArrayBoxFactory const*>{fac});
            linop->setSigma(sigma);
            linop->setDomainBC(lobc, hibc);
            set_eb_dirichlet(*linop);
            auto mlmg = std::make_unique<amrex::MLMG>(*linop);
            applyMLMGOptions(*mlmg, options);
            mlmg->setConvergenceNormType(amrex::MLMGNormType::greater);
            // In the collocated case, computeE needs ghost node data: ask
            // MLMG to fill the BC after it solves the problem.
            mlmg->setFinalFillBC(ng > 0);
            mlmg->setBottomSolver(parms.bottom_solver);
            mlmg->setNoGpuSync(parms.no_gpu_sync);
            return Pair{std::move(linop), std::move(mlmg)};
        };

        if (rho && options.tune_mlmg) {
            tune(lev, grids, dmap, factory, info, options, make, *rho,
                 relative_tolerance, absolute_tolerance);
        }

        if (!m_mlmg) {
            // Not tuned, or no trial configuration converged
            auto tmp = make(grids, dmap, &factory, m_parms);
            m_linop = std::move(tmp.first);
            m_mlmg = std::move(tmp.second);
        }
    }

    amrex::MLMG& getMLMG () { return *m_mlmg; }

    amrex::MLEBNodeFDLaplacian& getLinOp () { return *m_linop; }

    Parms const& getParms () const { return m_parms; }

    //! Tuning was skipped because the trial solve needed no iteration; retry later.
    [[nodiscard]] bool tuningSkipped () const { return m_tuning_skipped; }

    //! Set the EB potential, a value or a callable, @see MLEBNodeFDLaplacian::setEBDirichlet
    template <typename F>
    void setEBDirichlet (F const& f) { m_linop->setEBDirichlet(f); }

    //! Apply the per-solve options (tolerances, verbosity, ...) without losing the tuned ones
    void applyOptions (MLMGOptions const& options)
    {
        applyMLMGOptions(*m_mlmg, options);
        m_mlmg->setBottomSolver(m_parms.bottom_solver);
        m_mlmg->setNoGpuSync(m_parms.no_gpu_sync);
    }

    void solve (amrex::MultiFab& phi, amrex::MultiFab const& rho,
                amrex::Real relative_tolerance, amrex::Real absolute_tolerance)
    {
        if (m_parms.single_box) {
            solve_single_box(m_geom, *m_mlmg, m_ng, phi, rho, m_phi_tmp, m_rho_tmp,
                             relative_tolerance, absolute_tolerance);
        } else {
            m_mlmg->solve({&phi}, {&rho}, relative_tolerance, absolute_tolerance);
        }
    }

private:

    using Pair = std::pair<std::unique_ptr<amrex::MLEBNodeFDLaplacian>,
                           std::unique_ptr<amrex::MLMG>>;

    Parms m_parms;
    amrex::Geometry m_geom;
    //! The single-box layout's EB factory and temporaries. MLMG aliases the
    //! temporaries, and MLMG::getGradSolution needs them after the solve.
    std::unique_ptr<amrex::EBFArrayBoxFactory> m_single_box_factory;
    amrex::MultiFab m_phi_tmp;
    amrex::MultiFab m_rho_tmp;
    std::unique_ptr<amrex::MLEBNodeFDLaplacian> m_linop;
    std::unique_ptr<amrex::MLMG> m_mlmg; //!< references m_linop, declared after it
    int m_ng = 0;
    bool m_tuning_skipped = false;

    //! Elapsed time marking a configuration that failed and must not be selected
    static constexpr double t_failed = std::numeric_limits<double>::max();

    static std::string describe (Parms const& p)
    {
        std::string msg;
        if (p.single_box) { msg += "single_box=true "; }
        msg.append("max_coarsening_level=").append(std::to_string(p.mcl));
        if (p.ags != -1) {
            msg.append(" agglomeration_grid_size=").append(std::to_string(p.ags));
        }
        msg.append(" bottom_solver=").append(bottomSolverName(p.bottom_solver));
#ifdef AMREX_USE_GPU
        msg.append(" no_gpu_sync=").append(p.no_gpu_sync ? "true" : "false");
#endif
        return msg;
    }

    /** Time trial solves over the parameter sets; keep the fastest in m_parms,
     *  and its operator and solver in m_linop and m_mlmg */
    template <typename MakeF>
    void tune (int lev,
               amrex::BoxArray const& grids,
               amrex::DistributionMapping const& dmap,
               amrex::EBFArrayBoxFactory const& factory,
               amrex::LPInfo const& info,
               MLMGOptions const& options,
               MakeF const& make,
               amrex::MultiFab const& rho,
               amrex::Real relative_tolerance,
               amrex::Real absolute_tolerance)
    {
        int const verbosity = options.verbosity;
        int const trial_verbosity = (verbosity > 2) ? verbosity - 2 : 0;
        Parms const base = m_parms;

        amrex::MultiFab phi(rho.boxArray(), rho.DistributionMap(), 1, 1,
                            amrex::MFInfo{}, rho.Factory());

        // Some of the configurations tried below might blow up.
        FPExceptGuard const fpe_guard{};

        Pair best;                // the fastest operator and solver so far
        double tsolve = t_failed; // and its time

        // Time one configuration: the fastest of a few solves after one
        // untimed warm-up solve, which builds the setup.  t_failed if a
        // solve fails or gives NaNs.  Stops early when not competitive.
        auto const time_config = [&] (auto&& solve_fn, double t_best) -> double
        {
            constexpr int n_timed = 3;
            double t = t_failed;
            try {
                for (int itry = 0; itry <= n_timed; ++itry) {
                    phi.setVal(0.0);
                    amrex::Gpu::synchronize();
                    auto const t0 = amrex::second();
                    solve_fn();
                    amrex::Gpu::synchronize();
                    double t_try = amrex::second() - t0;
                    if (phi.contains_nan(0, 1, 0)) {
                        throw amrex::MLMG::error("MLMGWrapper: solution contains NaNs");
                    }
                    amrex::ParallelAllReduce::Max(t_try, amrex::ParallelContext::CommunicatorSub());
                    if (itry == 0) { continue; } // warm-up
                    t = std::min(t, t_try);
                    if (t_best < t_failed && t_try > 2.0*t_best) { break; }
                }
            } catch (amrex::MLMG::error const&) {
                t = t_failed;
            }
            return t;
        };

        auto const prepare_trial = [&] (Pair& tmp, int max_iters) {
            tmp.second->setThrowException(true); // trial solve: failure is caught, not fatal
            tmp.second->setVerbose(trial_verbosity);
            tmp.second->setMaxIter(max_iters);
        };

        auto const consider = [&] (Parms const& p, double t) {
            if (t < tsolve) {
                tsolve = t;
                m_parms = p;
            }
            if (verbosity > 1) {
                amrex::Print() << "MLMG tuning: " << describe(p) << " : "
                               << ((t < t_failed) ? std::to_string(t) + " s" : "failed") << "\n";
            }
        };

        // Baseline: the parameters as given
        int mcl0 = 0;
        int niter0 = options.max_iters;
        {
            Pair tmp = make(grids, dmap, &factory, base);
            prepare_trial(tmp, options.max_iters);
            consider(base, time_config([&] {
                tmp.second->solve({&phi}, {&rho}, relative_tolerance, absolute_tolerance);
            }, t_failed));
            // The MG hierarchy is built at the first solve.
            mcl0 = tmp.first->NMGLevels(0) - 1;
            if (tsolve < t_failed) {
                niter0 = tmp.second->getNumIters();
                best = std::move(tmp);
            }
        }
        if (tsolve < t_failed && niter0 == 0) {
            // Nothing to time with this right-hand side, e.g. rho = 0 at initialization.
            m_tuning_skipped = true;
            adopt(best);
            return;
        }
        // A configuration needing far more iterations than the baseline is no candidate.
        int const trial_max_iters = std::min(options.max_iters, 2*niter0);

        // On GPUs, MLMG can run in a single-stream no-sync region, which
        // avoids stream synchronizations at the cost of the concurrency of
        // multiple streams.  On CPUs this makes no difference.
#ifdef AMREX_USE_GPU
        constexpr std::array<bool,2> sync_options{true, false};
#else
        constexpr std::array<bool,1> sync_options{true};
#endif
        auto const push_unique = [] (auto& v, auto const& x) {
            if (std::find(v.begin(), v.end(), x) == v.end()) { v.push_back(x); }
        };
        std::vector<amrex::BottomSolver> bottom_solvers{base.bottom_solver};
        for (auto const b : {amrex::BottomSolver::Default, amrex::BottomSolver::algmg,
                             amrex::BottomSolver::smoother}) {
            push_unique(bottom_solvers, b);
        }
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        // The operator's own bottom solver: a fused kernel when the bottom
        // level is one box, otherwise the same as bicgstab.
        push_unique(bottom_solvers, amrex::BottomSolver::custom);
#endif
        std::vector<int> agg_grid_sizes{base.ags};
        if (grids.size() > 1 && info.do_agglomeration) {
            int const d = amrex::LPInfo::getDefaultAgglomerationGridSize();
            for (int const a : {d, 2*d, 4*d}) { push_unique(agg_grid_sizes, a); }
        }
        // The smoother alone cannot solve a large bottom problem.
        auto const skip = [&] (amrex::BottomSolver b, int mcl) {
            return (b == amrex::BottomSolver::smoother) && (mcl != mcl0);
        };
        // Keep a trial pair if one of its configurations became the fastest
        auto const keep_if_best = [&] (Pair& tmp, double t_before) {
            if (tsolve < t_before) { best = std::move(tmp); }
        };

        // One box covering the whole domain.  Its EB data, temporaries and
        // operator all live on one rank, so only small domains are tried.
        if ((lev == 0) && (grids.size() > 1)
            && (m_geom.Domain().numPts() <= (amrex::Long(128)*amrex::Long(128)*amrex::Long(128))))
        {
            amrex::BoxArray const sb_ba(m_geom.Domain());
            amrex::DistributionMapping const sb_dm(amrex::Vector<int>{0});
            m_single_box_factory = amrex::makeEBFabFactory(factory.getEBIndexSpace(), m_geom,
                                                           sb_ba, sb_dm, {1,1,1},
                                                           amrex::EBSupport::full);
            m_phi_tmp.define(amrex::convert(sb_ba, amrex::IntVect(1)), sb_dm, 1, 1);
            m_phi_tmp.setVal(0.0); // ghost nodes at Dirichlet faces are never written by MLMG
            m_rho_tmp.define(amrex::convert(sb_ba, amrex::IntVect(1)), sb_dm, 1, 0);

            // With a single box, the agglomeration grid size does not matter.
            for (int mcl = 0; mcl <= mcl0; ++mcl) {
                Parms p = base;
                p.single_box = true;
                p.mcl = mcl;
                Pair tmp = make(sb_ba, sb_dm, m_single_box_factory.get(), p);
                prepare_trial(tmp, trial_max_iters);
                double const t_before = tsolve;
                for (auto const b : bottom_solvers) {
                    if (skip(b, mcl)) { continue; }
                    p.bottom_solver = b;
                    tmp.second->setBottomSolver(b);
                    for (bool const nosync : sync_options) {
                        p.no_gpu_sync = nosync;
                        tmp.second->setNoGpuSync(nosync);
                        consider(p, time_config([&] {
                            solve_single_box(m_geom, *tmp.second, m_ng, phi, rho,
                                             m_phi_tmp, m_rho_tmp,
                                             relative_tolerance, absolute_tolerance);
                        }, tsolve));
                    }
                }
                keep_if_best(tmp, t_before);
            }
        }

        // The original grids
        for (int mcl = 0; mcl <= mcl0; ++mcl) {
            for (int const ags : agg_grid_sizes) {
                // Without coarse levels there is no agglomeration.
                if ((mcl == 0) && (ags != base.ags)) { continue; }
                Parms p = base;
                p.mcl = mcl;
                p.ags = ags;
                Pair tmp = make(grids, dmap, &factory, p);
                prepare_trial(tmp, trial_max_iters);
                double const t_before = tsolve;
                for (auto const b : bottom_solvers) {
                    if (skip(b, mcl)) { continue; }
                    p.bottom_solver = b;
                    tmp.second->setBottomSolver(b);
                    for (bool const nosync : sync_options) {
                        p.no_gpu_sync = nosync;
                        // The baseline, already timed
                        if (mcl == mcl0 && ags == base.ags && b == base.bottom_solver
                            && nosync == base.no_gpu_sync) { continue; }
                        tmp.second->setNoGpuSync(nosync);
                        consider(p, time_config([&] {
                            tmp.second->solve({&phi}, {&rho},
                                              relative_tolerance, absolute_tolerance);
                        }, tsolve));
                    }
                }
                keep_if_best(tmp, t_before);
            }
        }

        if (tsolve < t_failed) {
            if (verbosity > 0) {
                amrex::Print() << "MLMG parameters chosen on level " << lev << ": "
                               << describe(m_parms) << "\n";
            }
        } else {
            // Nothing worked with this right-hand side. Fall back to the
            // given parameters and hope the solver works with the real inputs.
            m_parms = base;
            if (verbosity > 0) {
                amrex::Print() << "MLMG tuning failed on level " << lev
                               << ". Using the default parameters.\n";
            }
        }
        if (!m_parms.single_box) {
            m_single_box_factory.reset();
            m_phi_tmp.clear();
            m_rho_tmp.clear();
        }
        adopt(best);
    }

    //! Make a trial pair the production operator and solver; it keeps its setup
    void adopt (Pair& best)
    {
        if (best.second) {
            best.second->setThrowException(false);
            m_linop = std::move(best.first);
            m_mlmg = std::move(best.second);
        }
    }

    //! RAII guard disabling floating point exception trapping.
    struct FPExceptGuard
    {
        FPExceptGuard ()
            : m_prev(amrex::OpenMP::get_max_threads(), amrex::FPExcept::none)
        {
            // Each thread has its own settings.
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            {
                m_prev[amrex::OpenMP::get_thread_num()]
                    = amrex::disableFPExcept(amrex::FPExcept::all);
            }
        }

        ~FPExceptGuard ()
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            {
                std::feclearexcept(FE_ALL_EXCEPT);
                amrex::setFPExcept(m_prev[amrex::OpenMP::get_thread_num()]);
            }
        }

        FPExceptGuard (FPExceptGuard const&) = delete;
        FPExceptGuard (FPExceptGuard &&) = delete;
        FPExceptGuard& operator= (FPExceptGuard const&) = delete;
        FPExceptGuard& operator= (FPExceptGuard &&) = delete;
    private:
        amrex::Vector<amrex::FPExcept> m_prev;
    };

    //! Solve on a single box covering the whole domain.
    static void solve_single_box (amrex::Geometry const& geom,
                                  amrex::MLMG& mlmg, int ng,
                                  amrex::MultiFab& phi, amrex::MultiFab const& rho,
                                  amrex::MultiFab& phi_tmp, amrex::MultiFab& rho_tmp,
                                  amrex::Real relative_tolerance,
                                  amrex::Real absolute_tolerance)
    {
        phi_tmp.ParallelCopy(phi, 0, 0, 1);
        rho_tmp.ParallelCopy(rho, 0, 0, 1);
        mlmg.solve({&phi_tmp}, {&rho_tmp}, relative_tolerance, absolute_tolerance);
        int const ngsrc = ng ? 1 : 0;
        int const ngdst = std::min(ngsrc, phi.nGrow());
        phi.ParallelCopy(phi_tmp, 0, 0, 1, ngsrc, ngdst, geom.periodicity());
    }
};

namespace details
{

/** The cached solver of one level and what it was built for */
struct EBPoissonSolverEntry
{
    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
    amrex::EBFArrayBoxFactory const* factory = nullptr;
    amrex::LPInfo info;
    MLMGWrapper::BCArray lobc{};
    MLMGWrapper::BCArray hibc{};
    MLMGWrapper::RealArray sigma{};
    bool unstable = false; //!< the coefficients changed on this level: no reuse, no tuning
    std::shared_ptr<MLMGWrapper> solver;
};

/** The per-level EB solver cache; freed by amrex::Finalize */
inline amrex::Vector<EBPoissonSolverEntry>& ebPoissonSolverCache ()
{
    static amrex::Vector<EBPoissonSolverEntry> cache;
    static bool finalize_registered = false;
    if (!finalize_registered) {
        finalize_registered = true;
        // Release the AMReX data before the arenas go away.  AMReX may be
        // initialized again in the same process (e.g., from Python).
        amrex::ExecOnFinalize([] () { cache.clear(); finalize_registered = false; });
    }
    return cache;
}

inline bool sameLPInfo (amrex::LPInfo const& a, amrex::LPInfo const& b)
{
    return a.do_agglomeration == b.do_agglomeration
        && a.do_consolidation == b.do_consolidation
        && a.do_semicoarsening == b.do_semicoarsening
        && a.agg_grid_size == b.agg_grid_size
        && a.con_grid_size == b.con_grid_size
        && a.con_ratio == b.con_ratio
        && a.con_strategy == b.con_strategy
        && a.has_metric_term == b.has_metric_term
        && a.max_coarsening_level == b.max_coarsening_level
        && a.max_semicoarsening_level == b.max_semicoarsening_level
        && a.semicoarsening_direction == b.semicoarsening_direction
        && a.hidden_direction == b.hidden_direction;
}

inline bool sameGeometry (amrex::Geometry const& a, amrex::Geometry const& b)
{
    if (a.Domain() != b.Domain() || a.isPeriodic() != b.isPeriodic()) { return false; }
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (a.ProbLo(idim) != b.ProbLo(idim) || a.ProbHi(idim) != b.ProbHi(idim)) {
            return false;
        }
    }
    return true;
}

/** Return the EB operator and solver of level `lev`, from the cache when possible
 *
 * The cached solver is reused while the level's geometry, grids, EB factory,
 * coarsening options, domain BCs and coefficients are unchanged.  A level
 * whose coefficients change (a beta that varies between solves) gets a
 * solver per call and no tuning, since its matrix is not stable.
 *
 * Parameters are those of the MLMGWrapper constructor, plus
 * \param[out] reused true if the solver comes from the cache (its EB potential is stale)
 */
template <typename F>
std::shared_ptr<MLMGWrapper>
getEBPoissonSolver (
    int lev,
    amrex::Geometry const& geom,
    amrex::BoxArray const& grids,
    amrex::DistributionMapping const& dmap,
    amrex::EBFArrayBoxFactory const& factory,
    amrex::LPInfo const& info,
    MLMGWrapper::BCArray const& lobc,
    MLMGWrapper::BCArray const& hibc,
    MLMGWrapper::RealArray const& sigma,
    MLMGOptions const& options,
    int ng,
    F const& set_eb_dirichlet,
    amrex::MultiFab const& rho,
    amrex::Real relative_tolerance,
    amrex::Real absolute_tolerance,
    bool& reused)
{
    reused = false;
    auto& cache = ebPoissonSolverCache();
    if (lev >= static_cast<int>(cache.size())) { cache.resize(lev+1); }
    auto& c = cache[lev];

    bool const same_layout = c.solver
        && sameGeometry(c.geom, geom)
        && c.grids == grids
        && c.dmap == dmap
        && c.factory == &factory
        && sameLPInfo(c.info, info)
        && c.lobc == lobc
        && c.hibc == hibc;

    if (same_layout && c.sigma == sigma) {
        bool const retry_tuning = options.tune_mlmg && c.solver->tuningSkipped();
        if (!retry_tuning) {
            reused = true;
            return c.solver;
        }
    } else if (same_layout) {
        c.unstable = true;
        c.solver.reset();
    }

    if (c.unstable) {
        return std::make_shared<MLMGWrapper>(lev, geom, grids, dmap, factory, info, lobc, hibc,
                                             sigma, options, ng, set_eb_dirichlet, nullptr,
                                             relative_tolerance, absolute_tolerance);
    }

    c.solver.reset(); // before building the next one, to limit the memory use
    c.solver = std::make_shared<MLMGWrapper>(lev, geom, grids, dmap, factory, info, lobc, hibc,
                                             sigma, options, ng, set_eb_dirichlet, &rho,
                                             relative_tolerance, absolute_tolerance);
    c.geom = geom;
    c.grids = grids;
    c.dmap = dmap;
    c.factory = &factory;
    c.info = info;
    c.lobc = lobc;
    c.hibc = hibc;
    c.sigma = sigma;
    return c.solver;
}

} // namespace details
} // namespace ablastr::fields

#endif // AMREX_USE_EB

#endif // ABLASTR_MLMG_WRAPPER_H_
