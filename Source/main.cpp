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
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXrocfftUtil.H"

#include <ablastr/parallelization/MPIInitHelpers.H>
#include <ablastr/utils/timer/Timer.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_Print.H>

// OLCFDEV-1655: Segfault during MPI_Init & in PMI_Allgather
// https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html#olcfdev-1655-occasional-seg-fault-during-mpi-init
#if defined(AMREX_USE_HIP)
#include <hip/hip_runtime.h>
#endif

#include <iostream>


int main(int argc, char* argv[])
{
// OLCFDEV-1655: Segfault during MPI_Init & in PMI_Allgather
// https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html#olcfdev-1655-occasional-seg-fault-during-mpi-init
#if defined(AMREX_USE_HIP) && defined(AMREX_USE_MPI)
    hipError_t hip_ok = hipInit(0);
    if (hip_ok != hipSuccess) {
        std::cerr << "hipInit failed with error code " << hip_ok << "! Aborting now.\n";
        return 1;
    }
#endif

    ablastr::parallelization::mpi_init(argc, argv);

    warpx::initialization::amrex_init(argc, argv);

    utils::rocfft::setup();

    {
        WARPX_PROFILE_VAR("main()", pmain);

        auto timer = ablastr::utils::timer::Timer{};
        timer.record_start_time();

        auto& warpx = WarpX::GetInstance();

        warpx.InitData();

        warpx.Evolve();

        //Print warning messages at the end of the simulation
        amrex::Print() <<
            ablastr::warn_manager::GetWMInstance().PrintGlobalWarnings("THE END");

        timer.record_stop_time();
        if (warpx.Verbose()) {
            amrex::Print() << "Total Time                     : "
                    << timer.get_global_duration() << '\n';
        }

        WARPX_PROFILE_VAR_STOP(pmain);

        WarpX::Finalize();
    }

    utils::rocfft::cleanup();

    amrex::Finalize();

    ablastr::parallelization::mpi_finalize ();
}
