/* Copyright 2020 Luca Fedeli, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_schwinger_process_wrapper_h_
#define WARPX_schwinger_process_wrapper_h_

#include "QedWrapperCommons.H"

//#define PXRMP_CORE_ONLY allows importing only the 'core functions' of the
//Schwinger process engine of the QED PICSAR library.
#define PXRMP_CORE_ONLY
#include <schwinger_pair_engine.hpp>
#include <AMReX_Random.H>

#include <cmath>

AMREX_GPU_DEVICE amrex::Real
getSchwingerProductionNumber (const amrex::Real dVdt,
                  const amrex::Real Ex, const amrex::Real Ey, const amrex::Real Ez,
                  const amrex::Real Bx, const amrex::Real By, const amrex::Real Bz,
                  const amrex::Real PoissonToGaussianThreshold = 100.)
{
    using namespace amrex;
    using namespace picsar::multi_physics;

    //If SI units are used, PICSAR QED does not use a reference length. 
    // However, a "lambda" parameter must still be passed to some functions. 
    const auto dummy_lambda = 1._rt;

    const auto expectedPairNumber =   dVdt *
        schwinger_pair_engine<Real,QedUtils::DummyStruct>::
        internal_compute_schwinger_pair_production_rate(
            Ex, Ey, Ez, Bx, By, Bz, dummy_lambda);

    if (expectedPairNumber <= PoissonToGaussianThreshold) {
        return RandomPoisson(expectedPairNumber);
    }
    else {
        const auto numpairs = RandomNormal(expectedPairNumber,sqrt(expectedPairNumber));
        return numpairs > 0._rt ? numpairs : 0._rt;
    }

}

#endif // WARPX_schwinger_process_wrapper_h_
