/* Copyright 2022 The ABLASTR Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: David Grote, Axel Huebl, Luca Fedeli
 */

#include "SafeCast.H"

#include "TextMsg.H"  // for ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE

#include <cmath>
#include <limits>


namespace
{
    template< typename int_type >
    AMREX_FORCE_INLINE
    int_type safeCastTo (
        const amrex::Real x,
        const std::string& real_name
    )
    {
        auto result = int_type(0);
        bool error_detected = false;
        std::string assert_msg;
        // (2.0*(numeric_limits<int>::max()/2+1)) converts numeric_limits<int>::max()+1 to a real ensuring accuracy to all digits
        // This accepts x = 2**31-1 but rejects 2**31.
        using namespace amrex::literals;
        constexpr int_type half_max_plus_one = std::numeric_limits<int_type>::max()/2+1;
        constexpr amrex::Real max_range = (2.0_rt*static_cast<amrex::Real>(half_max_plus_one));
        if (x < max_range) {
            if (std::ceil(x) >= std::numeric_limits<int_type>::min()) {
                result = static_cast<int_type>(x);
            } else {
                error_detected = true;
                assert_msg = "Negative overflow detected when casting " + real_name + " = " +
                             std::to_string(x) + " to integer type";
            }
        } else if (x > 0) {
            error_detected = true;
            assert_msg =  "Overflow detected when casting " + real_name + " = " + std::to_string(x) + " to integer type";
        } else {
            error_detected = true;
            assert_msg =  "NaN detected when casting " + real_name + " to integer type";
        }
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(!error_detected, assert_msg);
        return result;
   }
}


int
ablastr::utils::safeCastToInt (
    const amrex::Real x,
    const std::string& real_name
)
{
    return ::safeCastTo<int> (x, real_name);
}


long
ablastr::utils::safeCastToLong (
    const amrex::Real x,
    const std::string& real_name
)
{
    return ::safeCastTo<long> (x, real_name);
}
