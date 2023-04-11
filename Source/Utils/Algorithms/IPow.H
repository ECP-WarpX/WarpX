/* Copyright 2023 Luca Fedeli
 * Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>

namespace utils::algorithms
{
    /**
    * \brief This function is a replacement for std::pow when an integer power
    * known at compile time is used.
    *
    * @tparam Power the integer exponent (must be >= 1)
    * @tparam T the type of the argument
    * @param[in] x the valule to be exponentiated
    */
    template<unsigned int IPower, typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    constexpr auto ipow(T x) noexcept
    {
        static_assert(IPower >= 1, "Power must be > 1!");

        if constexpr (IPower == 1)
            return x;
        else if constexpr (IPower == 2)
            return x*x;
        else if constexpr (IPower%2 == 0)
            return ipow<2>(ipow<IPower/2>(x));
        else
            return x * ipow <IPower-1>(x);
    }
}
