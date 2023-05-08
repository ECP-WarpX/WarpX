/* Copyright 2023 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpXrocfftUtil.H"

#include <AMReX_Config.H>

#if defined(AMREX_USE_HIP) && defined(WARPX_USE_PSATD)
#  if __has_include(<rocfft/rocfft.h>)  // ROCm 5.3+
#    include <rocfft/rocfft.h>
#  else
#    include <rocfft.h>
#  endif
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
