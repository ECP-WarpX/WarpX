/* Copyright 2020 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MATH_UTILS_
#define WARPX_MATH_UTILS_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

namespace utils{

    /* \brief calculates the cross product between two vectors (must have size >=3)
    *
    * \tparam VectorType: the vector type
    * \param va: the first vector
    * \param vb: the second vector
    * \return a VectorType v1 x v2
    */
    template<typename VectorType>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    VectorType CrossProduct (const VectorType& va, const VectorType& vb)
    {
        return {
            va[1]*vb[2]-va[2]*vb[1],
            va[2]*vb[0]-va[0]*vb[2],
            va[0]*vb[1]-va[1]*vb[0] };
    }

}

#endif //WARPX_MATH_UTILS_
