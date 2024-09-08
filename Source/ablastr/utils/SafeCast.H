/* Copyright 2022 The ABLASTR Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: David Grote, Axel Huebl, Luca Fedeli
 */

#ifndef ABLASTR_UTILS_SAFE_CAST_H_
#define ABLASTR_UTILS_SAFE_CAST_H_

#include <AMReX_REAL.H>

#include <string>


namespace ablastr::utils
{
    /**
    * \brief Do a safe cast of a real to an int
    * This ensures that the float value is within the range of ints and if not,
    * raises an exception.
    *
    * \param x Real value to cast
    * \param real_name String, the name of the variable being casted to use in the error message
    */
    int
    safeCastToInt (amrex::Real x, const std::string& real_name);


    /**
    * \brief Do a safe cast of a real to a long
    * This ensures that the float value is within the range of longs and if not,
    * raises an exception.
    *
    * \param x Real value to cast
    * \param real_name String, the name of the variable being casted to use in the error message
    */
    long
    safeCastToLong (amrex::Real x, const std::string& real_name);
}

#endif // ABLASTR_UTILS_SAFE_CAST_H_
