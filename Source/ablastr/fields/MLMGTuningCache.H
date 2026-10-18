/* Copyright 2026 The WarpX Community
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_MLMG_TUNING_CACHE_H_
#define ABLASTR_MLMG_TUNING_CACHE_H_

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_Vector.H>

#include <map>
#include <memory>

namespace ablastr::fields
{

/** Settings and shared state of the MLMG parameter tuning, @see MLMGWrapper
 *
 * The maps are cleared by amrex::ExecOnFinalize.
 */
struct MLMGTuningCache
{
    /** The MLMG parameters that the tuning chooses from */
    struct Parms {
        bool single_box = false; //!< solve on one box covering the whole domain
        int mcl = 30;            //!< maximum coarsening level (AMReX default)
        int ags = -1;            //!< agglomeration grid size, -1 for the AMReX default
        bool no_gpu_sync = true; //!< GPU only: single stream, no stream synchronization
        amrex::BottomSolver bottom_solver = amrex::BottomSolver::Default;
    };

    //! Enable/disable tuning (opt-in, since it performs trial solves).
    static void setTuning (bool enable) { m_tune = enable; }
    static bool tuning () { return m_tune; }

    //! Layout of the single-box solve, plus its EB factory (built by the caller on demand).
    struct SingleBox {
        amrex::Box domain;
        amrex::BoxArray ba; //!< cell-centered
        amrex::DistributionMapping dm;
        std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> eb_factory;
    };

    //! Single-box layout of a level. Rebuilt if the domain changed.
    static SingleBox& singleBox (int lev, amrex::Box const& domain)
    {
        registerFinalize();
        auto& sb = m_single_box_map[lev];
        if (sb.ba.empty() || sb.domain != domain) {
            sb = SingleBox{};
            sb.domain = domain;
            sb.ba = amrex::BoxArray(domain);
            sb.dm = amrex::DistributionMapping(amrex::Vector<int>{0});
        }
        return sb;
    }

    //! Drop a level's single-box layout and EB factory, e.g. when single-box was not chosen.
    static void eraseSingleBox (int lev) { m_single_box_map.erase(lev); }

private:
    //! Register the clean-up once per amrex::Initialize.
    static void registerFinalize ()
    {
        if (!m_finalize_registered) {
            m_finalize_registered = true;
            amrex::ExecOnFinalize([] () {
                m_single_box_map.clear();
                m_finalize_registered = false;
            });
        }
    }

    inline static bool m_tune = false;
    inline static bool m_finalize_registered = false;
    inline static std::map<int, SingleBox> m_single_box_map;
};

} // namespace ablastr::fields

#endif // ABLASTR_MLMG_TUNING_CACHE_H_
