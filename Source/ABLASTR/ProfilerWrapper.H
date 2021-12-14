/* Copyright 2020-2021 Axel Huebl, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef ABLASTR_PROFILERWRAPPER_H_
#define ABLASTR_PROFILERWRAPPER_H_

#include <AMReX_BLProfiler.H>
#include <AMReX_GpuDevice.H>


namespace ablastr {

    AMREX_FORCE_INLINE
    void doDeviceSynchronize(bool const do_device_synchronize = false) {
        if (do_device_synchronize)
            amrex::Gpu::synchronize();
    }

    // Note that objects are destructed in the reverse order of declaration
    struct SynchronizeOnDestruct {
        SynchronizeOnDestruct(bool const do_device_synchronize = false)
                : m_do_device_synchronize(do_device_synchronize) {}

        AMREX_FORCE_INLINE
        ~SynchronizeOnDestruct() {
            doDeviceSynchronize(m_do_device_synchronize);
        }

        bool m_do_device_synchronize = false;
    };

} // namespace ablastr

// `BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__)` and `SYNC_V_##vname` used to make unique names for
// synchronizeOnDestruct objects, like `SYNC_SCOPE_0` and `SYNC_V_pmain`
#define ABLASTR_PROFILE(fname, sync) ablastr::doDeviceSynchronize(sync); BL_PROFILE(fname); ablastr::SynchronizeOnDestruct BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__){sync}
#define ABLASTR_PROFILE_VAR(fname, vname, sync) ablastr::doDeviceSynchronize(sync); BL_PROFILE_VAR(fname, vname); ablastr::SynchronizeOnDestruct SYNC_V_##vname{sync}
#define ABLASTR_PROFILE_VAR_NS(fname, vname, sync) BL_PROFILE_VAR_NS(fname, vname); ablastr::SynchronizeOnDestruct SYNC_V_##vname{sync}
#define ABLASTR_PROFILE_VAR_START(vname, sync) ablastr::doDeviceSynchronize(sync); BL_PROFILE_VAR_START(vname)
#define ABLASTR_PROFILE_VAR_STOP(vname, sync) ablastr::doDeviceSynchronize(sync); BL_PROFILE_VAR_STOP(vname)
#define ABLASTR_PROFILE_REGION(rname, sync) ablastr::doDeviceSynchronize(sync); BL_PROFILE_REGION(rname); ablastr::SynchronizeOnDestruct BL_PROFILE_PASTE(SYNC_R_, __COUNTER__){sync}

#endif // ABLASTR_PROFILERWRAPPER_H_
