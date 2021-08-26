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
#include <AMReX_Config.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_TinyProfiler.H>
#include <AMReX_Utility.H>

#if defined(AMREX_USE_MPI)
#  include <mpi.h>
#endif

#if defined(AMREX_USE_HIP) && defined(WARPX_USE_PSATD)
// cstddef: work-around for ROCm/rocFFT <=4.3.0
// https://github.com/ROCmSoftwarePlatform/rocFFT/blob/rocm-4.3.0/library/include/rocfft.h#L36-L42
#  include <cstddef>
#  include <rocfft.h>
#endif

int main(int argc, char* argv[])
{
    using namespace amrex;

    auto mpi_thread_levels = utils::warpx_mpi_init(argc, argv);

    warpx_amrex_init(argc, argv);

    utils::warpx_check_mpi_thread_level(mpi_thread_levels);

#if defined(AMREX_USE_HIP) && defined(WARPX_USE_PSATD)
    rocfft_setup();
#endif

    ParseGeometryInput();

    ConvertLabParamsToBoost();
    ReadBCParams();

#ifdef WARPX_DIM_RZ
    CheckGriddingForRZSpectral();
#endif

    WARPX_PROFILE_VAR("main()", pmain);

    const auto strt_total = static_cast<Real>(amrex::second());

    {
        WarpX warpx;

        warpx.RecordWarning("Topic1", "test_msg");
        warpx.RecordWarning("Topic1", "test_msg");
        warpx.RecordWarning("Topic1", "test_msg_2");
        warpx.RecordWarning("Topic2", "test_msg_2");
        warpx.RecordWarning("Topic2", "test_msg_2", WarnPriority::high);
        warpx.RecordWarning("Topic2", "test_msg_2", WarnPriority::low);
        warpx.RecordWarning("Topic3", "test_msg_3", WarnPriority::low);
        warpx.RecordWarning("Topic3", "test_msg_3", WarnPriority::low);
        warpx.RecordWarning("Topic3", "test_msg_3", WarnPriority::low);
        warpx.RecordWarning("MultiLine", "A\nmultiline\n\nmessage\ntest test test\n test test", WarnPriority::medium);
        warpx.RecordWarning("MultiLine2", "A\nmultiline\n\nmessage\ntest test test\n test test\n"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "a a a a a a a a a a a a a a a b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b b"
        "c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c \n c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c ", WarnPriority::medium);

        warpx.RecordWarning("Multiline", "first\nsecond\nthird\n\ndouble_space!!", WarnPriority::low);

        warpx.InitData();

        warpx.Evolve();

        warpx.PrintLocalWarnings("THE VERY END");

        if (warpx.Verbose()) {
            auto end_total = static_cast<Real>(amrex::second()) - strt_total;
            ParallelDescriptor::ReduceRealMax(end_total, ParallelDescriptor::IOProcessorNumber());
            Print() << "Total Time                     : " << end_total << '\n';
        }
    }

    WARPX_PROFILE_VAR_STOP(pmain);

#if defined(AMREX_USE_HIP) && defined(WARPX_USE_PSATD)
    rocfft_cleanup();
#endif

    Finalize();
#if defined(AMREX_USE_MPI)
    MPI_Finalize();
#endif
}
