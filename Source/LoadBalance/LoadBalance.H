/* Copyright 2024 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Luca Fedeli, Maxence Thevenet, Michael Rowan
 * Remi Lehe, Weiqun Zhang, levinem, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_LOADBALANCE_LOADBALANCE_H_
#define WARPX_LOADBALANCE_LOADBALANCE_H_

#include <AMReX_LayoutData.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <cstddef>

namespace warpx::load_balance
{
    /** Strategy to compute weights for use in load balance.
    */
    enum class CostsUpdateAlgo {
        Timers   , //!< load balance according to in-code timer-based weights (i.e., with  `costs`)
        Heuristic  /**< load balance according to weights computed from number of cells
                             and number of particles per box (i.e., with `costs_heuristic`) */
    };

    class CostTracker
    {
    public:
        CostTracker (int lev, std::size_t mfi_iter_index);
        ~CostTracker ();

    private:
        int m_lev;
        std::size_t m_mfi_iter_index;
        amrex::Real m_wt;
    };

    class AllCosts
    {
        friend class CostTracker;

    public:
        static AllCosts& get_instance()
        {
            static auto instance = AllCosts{};
            return instance;
        }

        AllCosts (AllCosts& other) = delete;
        void operator= (const AllCosts& other) = delete;

        [[nodiscard]]
        CostsUpdateAlgo get_update_algo () const noexcept
        {
            return m_update_algo;
        }

    protected:
        AllCosts ();
        ~ AllCosts (){}

    private:
        CostsUpdateAlgo m_update_algo;
        amrex::Vector<std::unique_ptr<amrex::LayoutData<amrex::Real> > > m_costs;
    };
}


#endif //WARPX_LOADBALANCE_LOADBALANCE_H_
