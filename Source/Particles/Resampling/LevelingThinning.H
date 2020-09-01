/* Copyright 2019-2020 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_LEVELING_THINNING_H_
#define WARPX_LEVELING_THINNING_H_

#include "Resampling.H"

/**
 * \brief This class will soon implement the leveling thinning algorithm, but for now it is just
 * a dummy class that prints "Resampling." used to test the resampling structure.
 */
class LevelingThinning: public ResamplingAlgorithm {
public:
    LevelingThinning ();

    ~LevelingThinning () override final;

    void operator() (WarpXParIter& pti, const int lev) const override final;

private:
    amrex::Real m_target_ratio;
};


#endif //WARPX_LEVELING_THINNING_H_
