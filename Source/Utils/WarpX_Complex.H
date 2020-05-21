/* Copyright 2019-2020 Andrew Myers, David Grote, Maxence Thevenet
 * Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMPLEX_H_
#define WARPX_COMPLEX_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuComplex.H>

#include <complex>

// Defines a complex type on GPU & CPU
using Complex = amrex::GpuComplex<amrex::Real>;

#ifdef AMREX_USE_GPU
#    ifdef WARPX_USE_PSATD
#    include <cufft.h>
#        ifdef AMREX_USE_FLOAT
static_assert( sizeof(Complex) == sizeof(cuComplex),
    "The complex types in WarpX and cuFFT do not match.");
#        else
static_assert( sizeof(Complex) == sizeof(cuDoubleComplex),
    "The complex types in WarpX and cuFFT do not match.");
#        endif
#    endif // WARPX_USE_PSATD

#else

#    ifdef WARPX_USE_PSATD
#    include <fftw3.h>
#        ifdef AMREX_USE_FLOAT
static_assert( sizeof(Complex) == sizeof(fftwf_complex),
    "The complex types in WarpX and FFTW do not match.");
#        else
static_assert( sizeof(Complex) == sizeof(fftw_complex),
    "The complex types in WarpX and FFTW do not match.");
#        endif
#    endif // WARPX_USE_PSATD

#endif // AMREX_USE_GPU
static_assert(sizeof(Complex) == sizeof(amrex::Real[2]),
              "Unexpected complex type.");

#endif //WARPX_COMPLEX_H_
