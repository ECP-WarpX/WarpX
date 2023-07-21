/* Copyright 2022 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_UTILS_ALGORITHMS_UPPER_BOUND_H_
#define WARPX_UTILS_ALGORITHMS_UPPER_BOUND_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

namespace utils::algorithms
{

    /** \brief Returns a pointer to the first element in the range [first, last) that is greater than val
     *
     * A re-implementation of the upper_bound algorithm suitable for GPU kernels.
     *
     * @param first: pointer to left limit of the range to consider
     * @param last: pointer to right limit of the range to consider
     * @param val: value to compare the elements of [first, last) to
     */
    template<typename T> AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    const T* upper_bound(const T* first, const T* last, const T& val)
    {
        const T* it;
        size_t count, step;
        count = last-first;
        while(count>0){
            it = first;
            step = count/2;
            it += step;
            if (!(val<*it)){
                first = ++it;
                count -= step + 1;
            }
            else{
                count = step;
            }
        }

        return first;
    }

}

#endif //WARPX_UTILS_ALGORITHMS_UPPER_BOUND_H_
