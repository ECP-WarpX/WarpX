/* Copyright 2020-2021 Axel Huebl, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PROFILERWRAPPER_H_
#define WARPX_PROFILERWRAPPER_H_

#include "ablastr/profiler/ProfilerWrapper.H"

namespace warpx::profiler
{
    /**
    * \brief
    * This singleton class holds the settings of the WarpX profiler
    * (for the moment the only setting is the flag storing the value of "do_device_synchronize" parameter)
    */
    class ProfileSettings
    {
        public:

        /**
        * \brief
        * This method returns a const references to the unique ProfileSettings instance.
        */
        static const ProfileSettings& GetInstance ();

        /**
        * \brief
        * This method initialises a new instance of the ProfileSettings.
        * The method throws an error if ProfileSettings is already initialized.
        */
        static void InitProfileSettings (bool device_synchronize_flag);

        /** Destructor */
        ~ProfileSettings () {delete m_instance;}

        /** Copy constructor */
        ProfileSettings ( ProfileSettings const &) = delete;
        /** Copy operator */
        ProfileSettings& operator= ( ProfileSettings const & ) = delete;
        /** Move constructor */
        ProfileSettings ( ProfileSettings && ) = default;
        /** Move operator */
        ProfileSettings& operator= ( ProfileSettings && ) = default;

        /**
        * \brief
        * This function returns the value of the flag "m_device_synchronize_flag"
        */
        bool get_device_synchronize_flag () const
        {
            return m_device_synchronize_flag;
        }

        private:

        /** Constructor */
        ProfileSettings (const bool device_synchronize_flag):
            m_device_synchronize_flag{device_synchronize_flag}
        {}

        /** Pointer to the (unique) ProfileSettings instance */
        static ProfileSettings* m_instance;

        bool m_device_synchronize_flag; /*whether to call a amrex::Gpu::synchronize() around profiling regions.*/

    };

    /**
    * \brief
    * Helper function to simplify the call to ProfileSettings::GetInstance().get_device_synchronize_flag()
    */
    bool get_device_synchronize_flag ();
}

// `BL_PROFILE_PASTE(SYNC_SCOPE_, __COUNTER__)` and `SYNC_V_##vname` used to make unique names for
// synchronizeOnDestruct objects, like `SYNC_SCOPE_0` and `SYNC_V_pmain`
#define WARPX_PROFILE(fname) ABLASTR_PROFILE(fname, warpx::profiler::get_device_synchronize_flag())
#define WARPX_PROFILE_VAR(fname, vname) ABLASTR_PROFILE_VAR(fname, vname, warpx::profiler::get_device_synchronize_flag())
#define WARPX_PROFILE_VAR_NS(fname, vname) ABLASTR_PROFILE_VAR_NS(fname, vname, warpx::profiler::get_device_synchronize_flag())
#define WARPX_PROFILE_VAR_START(vname) ABLASTR_PROFILE_VAR_START(vname, warpx::profiler::get_device_synchronize_flag())
#define WARPX_PROFILE_VAR_STOP(vname) ABLASTR_PROFILE_VAR_STOP(vname, warpx::profiler::get_device_synchronize_flag())
#define WARPX_PROFILE_REGION(rname) ABLASTR_PROFILE_REGION(rname, warpx::profiler::get_device_synchronize_flag())

#endif // WARPX_PROFILERWRAPPER_H_
