/* Copyright 2020-2021 Axel Huebl, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PROFILERWRAPPER_H_
#define WARPX_PROFILERWRAPPER_H_

#include "WarpX.H"
#include "ablastr/profiler/ProfilerWrapper.H"


// `BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__)` and `SYNC_V_##vname` used to make unique names for
// synchronizeOnDestruct objects, like `SYNC_SCOPE_0` and `SYNC_V_pmain`
#define WARPX_PROFILE(fname) ABLASTR_PROFILE(fname, WarpX::do_device_synchronize)
#define WARPX_PROFILE_VAR(fname, vname) ABLASTR_PROFILE_VAR(fname, vname, WarpX::do_device_synchronize)
#define WARPX_PROFILE_VAR_NS(fname, vname) ABLASTR_PROFILE_VAR_NS(fname, vname, WarpX::do_device_synchronize)
#define WARPX_PROFILE_VAR_START(vname) ABLASTR_PROFILE_VAR_START(vname, WarpX::do_device_synchronize)
#define WARPX_PROFILE_VAR_STOP(vname) ABLASTR_PROFILE_VAR_STOP(vname, WarpX::do_device_synchronize)
#define WARPX_PROFILE_REGION(rname) ABLASTR_PROFILE_REGION(rname, WarpX::do_device_synchronize)

#endif // WARPX_PROFILERWRAPPER_H_
