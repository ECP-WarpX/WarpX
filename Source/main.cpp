/* Copyright 2016-2020 Andrew Myers, Ann Almgren, Axel Huebl
 *                     David Grote, Jean-Luc Vay, Remi Lehe
 *                     Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Initialization/WarpXAMReXInit.H"
#include "Utils/MPIInitHelpers.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXrocfftUtil.H"
#include "Utils/WarpXUtil.H"

#include <ablastr/warn_manager/WarnManager.H>
#include <ablastr/utils/timer/Timer.H>

#include <AMReX_Print.H>

int main(int argc, char* argv[])
{
    utils::warpx_mpi_init(argc, argv);

    warpx_amrex_init(argc, argv);

    utils::rocfft::setup();

    ParseGeometryInput();

    ConvertLabParamsToBoost();
    ReadBCParams();

#ifdef WARPX_DIM_RZ
    CheckGriddingForRZSpectral();
#endif

    {
        WARPX_PROFILE_VAR("main()", pmain);

        auto timer = ablastr::utils::timer::Timer{};
        timer.record_start_time();

        WarpX warpx;

        warpx.InitData();

        warpx.Evolve();

        //Print warning messages at the end of the simulation
        ablastr::warn_manager::GetWMInstance().PrintGlobalWarnings("THE END");

        timer.record_stop_time();
        if (warpx.Verbose()) {
            amrex::Print() << "Total Time                     : "
                    << timer.get_global_duration() << '\n';
        }

        WARPX_PROFILE_VAR_STOP(pmain);
    }

    utils::rocfft::cleanup();

    amrex::Finalize();

    utils::warpx_mpi_finalize ();
}
