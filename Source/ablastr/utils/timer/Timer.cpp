/* Copyright 2023 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Timer.H"

#include <AMReX_ParallelDescriptor.H>

using namespace ablastr::utils::timer;

Timer::Timer(){}

void
Timer::record_start_time() noexcept
{
    m_start_time = amrex::ParallelDescriptor::second();
}

void
Timer::record_stop_time() noexcept
{
    m_stop_time = amrex::ParallelDescriptor::second();
}

amrex::Real
Timer::get_duration () noexcept
{
    return m_stop_time - m_start_time;
}

amrex::Real
Timer::get_global_duration ()
{
    auto duration = this->get_duration();
    amrex::ParallelDescriptor::ReduceRealMax(
        duration,
        amrex::ParallelDescriptor::IOProcessorNumber());
    return duration;
}
