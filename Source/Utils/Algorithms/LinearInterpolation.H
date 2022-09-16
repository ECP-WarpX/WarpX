/* Copyright 2022 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_UTILS_ALGORITHMS_LINEAR_INTERPOLATION_H_
#define WARPX_UTILS_ALGORITHMS_LINEAR_INTERPOLATION_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

namespace utils::algorithms
{
    /** \brief Performs a linear interpolation
     *
     * Performs a linear interpolation at x given the 2 points
     * (x0, f0) and (x1, f1)
     */
    template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T linear_interp(T x0, T x1, T f0, T f1, T x)
    {
        return ((x1-x)*f0 + (x-x0)*f1)/(x1-x0);
    }

    /** \brief Performs a bilinear interpolation
     *
     * Performs a bilinear interpolation at (x,y) given the 4 points
     * (x0, y0, f00), (x0, y1, f01), (x1, y0, f10), (x1, y1, f11).
     */
    template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T bilinear_interp(T x0, T x1, T y0, T y1, T f00, T f01, T f10, T f11, T x, T y)
    {
        const T fx0 = linear_interp(x0, x1, f00, f10, x);
        const T fx1 = linear_interp(x0, x1, f01, f11, x);
        return linear_interp(y0, y1, fx0, fx1, y);
    }

    /** \brief Performs a trilinear interpolation
     *
     * Performs a trilinear interpolation at (x,y,z) given the 8 points
     * (x0, y0, z0, f000), (x0, y0, z1, f001), (x0, y1, z0, f010), (x0, y1, z1, f011),
     * (x1, y0, z0, f100), (x1, y0, z1, f101), (x1, y1, z0, f110), (x1, y1, z1, f111)
     */
    template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T trilinear_interp(T x0, T x1,T y0, T y1, T z0, T z1,
        T f000, T f001, T f010, T f011, T f100, T f101, T f110, T f111,
        T x, T y, T z)
    {
        const T fxy0 = bilinear_interp(
            x0, x1, y0, y1, f000, f010, f100, f110, x, y);
        const T fxy1 = bilinear_interp(
            x0, x1, y0, y1, f001, f011, f101, f111, x, y);
        return linear_interp(z0, z1, fxy0, fxy1, z);
    }
}

#endif //WARPX_UTILS_ALGORITHMS_LINEAR_INTERPOLATION_H_
