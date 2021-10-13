/* Copyright 2020 Axel Huebl, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PROFILERWRAPPER_H_
#define WARPX_PROFILERWRAPPER_H_

#include <AMReX_BLProfiler.H>
#include <AMReX_GpuDevice.H>
#include <WarpX.H>

template<int detail_level>
AMREX_FORCE_INLINE
void doDeviceSynchronize ()
{
    if ( WarpX::do_device_synchronize >= detail_level )
        amrex::Gpu::synchronize();
}

// Note that objects are destructed in the reverse order of declaration
template<int detail_level>
struct synchronizeOnDestruct {
    AMREX_FORCE_INLINE
    ~synchronizeOnDestruct () {
        doDeviceSynchronize<detail_level>();
    }
};

// `BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__)` and `SYNC_V_##vname` used to make unique names for
// synchronizeOnDestruct objects, like `SYNC_SCOPE_0` and `SYNC_V_pmain`
#define WARPX_PROFILE(fname) doDeviceSynchronize<1>(); BL_PROFILE(fname); synchronizeOnDestruct<1> BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__){}
#define WARPX_PROFILE_VAR(fname, vname) doDeviceSynchronize<1>(); BL_PROFILE_VAR(fname, vname); synchronizeOnDestruct<1> SYNC_V_##vname{}
#define WARPX_PROFILE_VAR_NS(fname, vname) BL_PROFILE_VAR_NS(fname, vname); synchronizeOnDestruct<1> SYNC_V_##vname{}
#define WARPX_PROFILE_VAR_START(vname) doDeviceSynchronize<1>(); BL_PROFILE_VAR_START(vname)
#define WARPX_PROFILE_VAR_STOP(vname) doDeviceSynchronize<1>(); BL_PROFILE_VAR_STOP(vname)
#define WARPX_PROFILE_REGION(rname) doDeviceSynchronize<1>(); BL_PROFILE_REGION(rname); synchronizeOnDestruct<1> BL_PROFILE_PASTE(SYNC_R_, __COUNTER__){}

#define WARPX_DETAIL_PROFILE(fname) doDeviceSynchronize<2>(); BL_PROFILE(fname); synchronizeOnDestruct<2> BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__){}
#define WARPX_DETAIL_PROFILE_VAR(fname, vname) doDeviceSynchronize<2>(); BL_PROFILE_VAR(fname, vname); synchronizeOnDestruct<2> SYNC_V_##vname{}
#define WARPX_DETAIL_PROFILE_VAR_NS(fname, vname) BL_PROFILE_VAR_NS(fname, vname); synchronizeOnDestruct<2> SYNC_V_##vname{}
#define WARPX_DETAIL_PROFILE_VAR_START(vname) doDeviceSynchronize<2>(); BL_PROFILE_VAR_START(vname)
#define WARPX_DETAIL_PROFILE_VAR_STOP(vname) doDeviceSynchronize<2>(); BL_PROFILE_VAR_STOP(vname)
#define WARPX_DETAIL_PROFILE_REGION(rname) doDeviceSynchronize<2>(); BL_PROFILE_REGION(rname); synchronizeOnDestruct<2> BL_PROFILE_PASTE(SYNC_R_, __COUNTER__){}

#endif // WARPX_PROFILERWRAPPER_H_
