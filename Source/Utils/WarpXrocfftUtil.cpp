/* Copyright 2023 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpXrocfftUtil.H"

#include <AMReX_Config.H>

#if defined(AMREX_USE_HIP) && defined(WARPX_USE_PSATD)
// cstddef: work-around for ROCm/rocFFT <=4.3.0
// https://github.com/ROCmSoftwarePlatform/rocFFT/blob/rocm-4.3.0/library/include/rocfft.h#L36-L42
#  include <cstddef>
#  include <rocfft.h>
#endif

void
utils::rocfft::setup()
{
#if defined(AMREX_USE_HIP) && defined(WARPX_USE_PSATD)
    rocfft_setup();
#endif
}

void
utils::rocfft::cleanup()
{
#if defined(AMREX_USE_HIP) && defined(WARPX_USE_PSATD)
    rocfft_cleanup();
#endif
}
