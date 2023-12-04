/* Copyright 2022 Thomas Clark
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_UTILS_MATH_COMPLEX__H_
#define WARPX_UTILS_MATH_COMPLEX__H_

#include <AMReX_Config.H>
#include <AMReX_REAL.H>

#if defined(AMREX_USE_CUDA)
#   include <cuComplex.h>
#elif defined(AMREX_USE_HIP)
#   include <hip/hip_complex.h>    

#else
#   if WARPX_USE_PSATD
#       include <fftw3.h>   
#    else 
#       include <complex>

#   endif 
#endif

    // First, define library-dependent types (complex, FFT plan)
namespace warpx::utils
{
    /** Complex type for FFT, depends on FFT library */
#if defined(AMREX_USE_CUDA)
#  ifdef AMREX_USE_FLOAT
    using Complex = cuComplex;
#  else
    using Complex = cuDoubleComplex;
#  endif
#elif defined(AMREX_USE_HIP)
#  ifdef AMREX_USE_FLOAT
    using Complex = float2;
#  else
    using Complex = double2;
#  endif
#else
#if     WARPX_USE_PSATD
#         ifdef AMREX_USE_FLOAT
            using Complex = fftwf_complex;
#         else
            using Complex = fftw_complex;
#         endif
#         else 
            using Complex = std::complex<amrex::Real>;
#       endif    
#endif
}
#endif // WARPX_UTILS_MATH_COMPLEX__H_