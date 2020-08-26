/* Copyright 2019-2020 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_RESAMPLING_TRIGGER_H_
#define WARPX_RESAMPLING_TRIGGER_H_

#include "Utils/IntervalsParser.H"

#include <AMReX_ParmParse.H>

class ResamplingTrigger
{
public:

    ResamplingTrigger ();

    bool triggered (const int n, const amrex::Real global_numparts);

    void initialize_global_numcells ();

private:
    IntervalsParser m_resampling_intervals;
    amrex::Real m_max_avg_ppc = std::numeric_limits<amrex::Real>::max();
    amrex::Real m_global_numcells;
    bool m_initialized;
};

#endif //WARPX_RESAMPLING_TRIGGER_H_
