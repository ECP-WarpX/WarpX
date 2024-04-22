/* Copyright 2020-2024 Axel Huebl, Maxence Thevenet, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpXProfilerWrapper.H"

#include "Utils/TextMsg.H"

#include "AMReX_Extension.H"

using namespace warpx::profiler;

ProfileSettings* ProfileSettings::m_instance = nullptr;

const ProfileSettings& ProfileSettings::GetInstance ()
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_instance != nullptr,
        "warpx::profiler::ProfileSettings is not initialized.");
    AMREX_ASSUME(m_instance != nullptr);
    return *m_instance;
}


void ProfileSettings::InitProfileSettings (const bool device_synchronize_flag)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_instance == nullptr,
        "warpx::profiler::ProfileSettings is already initialized.");
    m_instance = new ProfileSettings{device_synchronize_flag};
}

bool warpx::profiler::get_device_synchronize_flag()
{
    return ProfileSettings::GetInstance().get_device_synchronize_flag();
}
