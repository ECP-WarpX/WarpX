/* Copyright 2024 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Luca Fedeli, Maxence Thevenet, Michael Rowan
 * Remi Lehe, Weiqun Zhang, levinem, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "LoadBalance.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_Gpu.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_ParmParse.H>

#include <string>

using namespace warpx::load_balance;
using namespace amrex::literals;

namespace
{
    CostsUpdateAlgo parse_cost_update_algo ()
    {
        using namespace std;

        const amrex::ParmParse pp_algo("algo");

        if (std::string buf; pp_algo.query("load_balance_costs_update", buf))
        {
            if (buf == "timers"s){
                return CostsUpdateAlgo::Timers;
            }
            else if (buf == "heuristic"s){
                return CostsUpdateAlgo::Heuristic;
            }
            else{
            WARPX_ABORT_WITH_MESSAGE(
                "'"s + buf + "'"s + " is not a valid cost update algorithm. "s +
                "Please select either 'heuristic' or 'timers'.");
            }
        }

        return CostsUpdateAlgo::Timers;
    }
}

CostTracker::CostTracker (int lev, std::size_t mfi_iter_index):
    m_lev{lev},
    m_mfi_iter_index{mfi_iter_index},
    m_wt{0.0_rt}
{
    if (AllCosts::get_instance().get_update_algo() == CostsUpdateAlgo::Heuristic)
        return;

    amrex::Gpu::synchronize();

    m_wt = static_cast<amrex::Real>(amrex::second());
}

void
CostTracker::add () const noexcept
{
    if (AllCosts::get_instance().get_update_algo() == CostsUpdateAlgo::Heuristic)
        return;

    amrex::Gpu::synchronize();
    const auto time_diff = static_cast<amrex::Real>(amrex::second()) - m_wt;
    auto& cost = AllCosts::get_instance().m_costs[m_lev];
    amrex::HostDevice::Atomic::Add( &(*cost)[m_mfi_iter_index], time_diff);
}

AllCosts::AllCosts ()
{
    m_update_algo = ::parse_cost_update_algo();
}
