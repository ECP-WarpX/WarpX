/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_MATH_UTILS_H_
#define ABLASTR_MATH_UTILS_H_

#include <AMReX_Config.H>
#include <AMReX_GpuQualifiers.H>

namespace ablastr::math
{
    template <typename T, typename Iter>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void LinSpaceFill(T start, T end, int size, Iter it)
    {
        const auto delta = static_cast<T>((end - start)/(size-1));

        for (int i = 0; i < size-1; ++i)
        {
            *(it++) = start + static_cast<T>(i*delta);
        }
        *it = end;
    }

}

#endif //ABLASTR_MATH_UTILS_H_
