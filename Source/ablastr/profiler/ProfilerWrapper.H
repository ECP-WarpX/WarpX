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


namespace ablastr::profiler
{
    /** Conditionally synchronizes active GPU operations
     *
     * @param do_device_synchronize perform amrex::Gpu::synchronize() if true
     */
    AMREX_FORCE_INLINE
    void
    device_synchronize(bool const do_device_synchronize = false) {
        if (do_device_synchronize) {
            amrex::Gpu::synchronize();
        }
    }

    /** An object that conditionally calls device_synchronize() on destruction
     *
     * Note that objects are destructed in the reverse order of declaration
     */
    struct SynchronizeOnDestruct {
        SynchronizeOnDestruct(bool const do_device_synchronize = false)
                : m_do_device_synchronize(do_device_synchronize) {}

        AMREX_FORCE_INLINE
        ~SynchronizeOnDestruct() {
            device_synchronize(m_do_device_synchronize);
        }

        // default move and copy operations
        SynchronizeOnDestruct(const SynchronizeOnDestruct&) = default;
        SynchronizeOnDestruct& operator=(const SynchronizeOnDestruct&) = default;
        SynchronizeOnDestruct(SynchronizeOnDestruct&&) = default;
        SynchronizeOnDestruct& operator=(SynchronizeOnDestruct&& field_data) = default;

        bool m_do_device_synchronize = false;
    };

} // namespace ablastr::profiler

// `BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__)` and `SYNC_V_##vname` used to make unique names for
// synchronizeOnDestruct objects, like `SYNC_SCOPE_0` and `SYNC_V_pmain`
#define ABLASTR_PROFILE(fname, sync) ablastr::profiler::device_synchronize(sync); BL_PROFILE(fname); const ablastr::profiler::SynchronizeOnDestruct BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__){sync}
#define ABLASTR_PROFILE_VAR(fname, vname, sync) ablastr::profiler::device_synchronize(sync); BL_PROFILE_VAR(fname, vname); const ablastr::profiler::SynchronizeOnDestruct SYNC_V_##vname{sync}
#define ABLASTR_PROFILE_VAR_NS(fname, vname, sync) BL_PROFILE_VAR_NS(fname, vname); const ablastr::profiler::SynchronizeOnDestruct SYNC_V_##vname{sync}
#define ABLASTR_PROFILE_VAR_START(vname, sync) ablastr::profiler::device_synchronize(sync); BL_PROFILE_VAR_START(vname)
#define ABLASTR_PROFILE_VAR_STOP(vname, sync) ablastr::profiler::device_synchronize(sync); BL_PROFILE_VAR_STOP(vname)
#define ABLASTR_PROFILE_REGION(rname, sync) ablastr::profiler::device_synchronize(sync); BL_PROFILE_REGION(rname); const ablastr::profiler::SynchronizeOnDestruct BL_PROFILE_PASTE(SYNC_R_, __COUNTER__){sync}

#endif // ABLASTR_PROFILERWRAPPER_H_
