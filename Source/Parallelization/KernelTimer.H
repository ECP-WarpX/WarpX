/* Copyright 2019-2020 Michael Rowan, Axel Huebl, Kevin Gott
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef KERNELTIMER_H_
#define KERNELTIMER_H_

#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_REAL.H>

// #ifdef WARPX_USE_GPUCLOCK

#include <climits>

/**
 * \brief Defines a timer object to be used on GPU; measures summed thread cycles.
 */
class KernelTimer
{
public:
    /** Constructor.
     * \param[in] do_timing Controls whether timer is active.
     * \param[in,out] cost Pointer to cost which holds summed thread cycles
     * (for performance, it is recommended to allocate pinned host memory).
     */
    AMREX_GPU_DEVICE
    KernelTimer
#if (defined AMREX_USE_GPU)
    (const bool do_timing, amrex::Real* cost)
        : m_do_timing(do_timing)
#else
    (const bool do_timing, amrex::Real* cost)
#endif
    {
#if (defined AMREX_USE_GPU)
#   if (defined WARPX_USE_GPUCLOCK)
        if (m_do_timing && cost) {
            m_cost = cost;
#       if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            // Start the timer
            m_wt = clock64();

#       elif defined(AMREX_USE_DPCPP)
            // To be updated
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( m_do_timing == false,
                                              "KernelTimer not yet supported for this hardware." );
#       endif
        }
#   else // WARPX_USE_GPUCLOCK
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( m_do_timing == false,
                                          "`algo.load_balance_costs_update = gpuclock` requires to compile with `-DWarpX_GPUCLOCK=ON`.");
#   endif // WARPX_USE_GPUCLOCK
#else  // AMREX_USE_GPU
        amrex::ignore_unused(do_timing, cost);
#endif // AMREX_USE_GPU
    }

    //! Destructor.
    AMREX_GPU_DEVICE
    ~KernelTimer () {
#if (defined AMREX_USE_GPU)
        if (m_do_timing && m_cost) {
#   if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            m_wt = clock64() - m_wt;
            amrex::Gpu::Atomic::Add( m_cost, amrex::Real(m_wt));
#   elif defined(AMREX_USE_DPCPP)
            // To be updated
#   endif
        }
#endif
    }

#if (defined AMREX_USE_GPU)
private:
    //! Stores whether kernel timer is active.
    bool m_do_timing;

    //! Location in which to accumulate costs from all threads.
    amrex::Real* m_cost = nullptr;

    //! Store the time difference (cost) from a single thread.
    long long int m_wt;
#endif
};
#endif
