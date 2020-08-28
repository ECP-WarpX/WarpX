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

/**
 * \brief This class will soon implement the leveling thinning algorithm, but for now it is just
 * a dummy class that prints "Resampling." used to test the resampling structure.
 */
class LevelingThinning: public ResamplingAlgorithm {
public:
    void operator() (WarpXParIter& /*pti*/) override {amrex::Print() << "Resampling. \n";}
};


#endif //WARPX_LEVELING_THINNING_H_
