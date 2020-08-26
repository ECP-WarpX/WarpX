/* Copyright 2019-2020 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_LEVELING_THINNING_H_
#define WARPX_LEVELING_THINNING_H_

#include "Resampling.H"

#include <AMReX_Print.H>

class LevelingThinning: public ResamplingAlgorithm {
public:
    void run(WarpXParIter& pti) override {amrex::Print() << "Resampling. \n";}
};


#endif //WARPX_LEVELING_THINNING_H_
