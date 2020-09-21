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
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>


int main(int argc, char* argv[])
{
    using namespace amrex;

    auto mpi_thread_levels = utils::warpx_mpi_init(argc, argv);

    warpx_amrex_init(argc, argv);

    utils::warpx_check_mpi_thread_level(mpi_thread_levels);

    ConvertLabParamsToBoost();

    CheckGriddingForRZSpectral();

    WARPX_PROFILE_VAR("main()", pmain);

    const Real strt_total = amrex::second();

    {
        WarpX warpx;

        warpx.InitData();

        warpx.Evolve();

        Real end_total = amrex::second() - strt_total;

        ParallelDescriptor::ReduceRealMax(end_total, ParallelDescriptor::IOProcessorNumber());
        if (warpx.Verbose()) {
            Print() << "Total Time                     : " << end_total << '\n';
            Print() << "WarpX Version: " << WarpX::Version() << '\n';
            Print() << "PICSAR Version: " << WarpX::PicsarVersion() << '\n';
        }
    }

    WARPX_PROFILE_VAR_STOP(pmain);

    Finalize();
#if defined(AMREX_USE_MPI)
    MPI_Finalize();
#endif
}
