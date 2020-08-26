/* Copyright 2019-2020 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_RESAMPLING_H_
#define WARPX_RESAMPLING_H_

#include "ResamplingTrigger.H"

class WarpXParIter;

struct ResamplingAlgorithm
{
    virtual void run(WarpXParIter& pti) {};
};

class Resampling
{
public:

    Resampling ();

    bool triggered (const int n, const amrex::Real global_numparts);

    void run (WarpXParIter& pti) const;

private:
    ResamplingTrigger m_resampling_trigger;
    std::unique_ptr<ResamplingAlgorithm> m_resampling_algorithm;
};

#endif //WARPX_RESAMPLING_H_
