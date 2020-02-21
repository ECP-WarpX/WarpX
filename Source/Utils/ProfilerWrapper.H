/* Copyright 2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PROFILER_WRAPPER_H_
#define WARPX_PROFILER_WRAPPER_H_

#include <AMReX_BLProfiler.H>

WX_PROFILE_VAR_NS WX_PROFILE_VAR_START WX_PROFILE_VAR_STOP WX_PROFILE  WX_PROFILE_REGION

void wrapper_profile_var_ns(std::string function_name, std::string variable_name,
                            int verbose=0, bool synchronize=false)
{
    WX_PROFILE_VAR_NS(function_name, blp_copy);
}

#endif // WARPX_PROFILER_WRAPPER_H_
