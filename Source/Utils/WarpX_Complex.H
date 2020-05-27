/* Copyright 2019-2020 Andrew Myers, David Grote, Maxence Thevenet
 * Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMPLEX_H_
#define WARPX_COMPLEX_H_

#include "FieldSolver/SpectralSolver/AnyFFT.H"

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuComplex.H>

#include <complex>

// Defines a complex type on GPU & CPU
using Complex = amrex::GpuComplex<amrex::Real>;

static_assert(sizeof(Complex) == sizeof(AnyFFT::Complex),
    "The complex type in WarpX and the FFT library do not match.");

static_assert(sizeof(Complex) == sizeof(amrex::Real[2]),
              "Unexpected complex type.");

#endif //WARPX_COMPLEX_H_
