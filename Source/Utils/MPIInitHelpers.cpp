/* Copyright 2020 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MPIInitHelpers.H"

#include "WarpX.H"

#include <AMReX_Config.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#if defined(AMREX_USE_MPI)
#   include <mpi.h>
#endif

#include <string>
#include <utility>
#include <sstream>

namespace utils
{
    int
    warpx_mpi_thread_required ()
    {
        int thread_required = -1;
#ifdef AMREX_USE_MPI
        thread_required = MPI_THREAD_SINGLE;  // equiv. to MPI_Init
#   ifdef AMREX_USE_OMP
        thread_required = MPI_THREAD_FUNNELED;
#   endif
#   ifdef AMREX_MPI_THREAD_MULTIPLE  // i.e. for async_io
        thread_required = MPI_THREAD_MULTIPLE;
#   endif
#endif
        return thread_required;
    }

    std::pair< int, int >
    warpx_mpi_init (int argc, char* argv[])
    {
        int thread_required = warpx_mpi_thread_required();
        int thread_provided = -1;
#ifdef AMREX_USE_MPI
        MPI_Init_thread(&argc, &argv, thread_required, &thread_provided);
#else
        amrex::ignore_unused(argc, argv);
#endif
        return std::make_pair(thread_required, thread_provided);
    }

    void
    warpx_check_mpi_thread_level ()
    {
#ifdef AMREX_USE_MPI
        int thread_required = warpx_mpi_thread_required();
        int thread_provided = -1;
        MPI_Query_thread(&thread_provided);
        auto mtn = amrex::ParallelDescriptor::mpi_level_to_string;

        std::stringstream ss;
        if( thread_provided < thread_required ){
            ss << "WARNING: Provided MPI thread safety level ("
                           << mtn(thread_provided) << ") is LOWER than requested "
                           << mtn(thread_required) << "). This might lead to undefined "
                           << "results in asynchronous operations (e.g. async_io).";
            WarpX::GetInstance().RecordWarning("MPI", ss.str(), WarnPriority::high);
        }
        if( thread_provided > thread_required ){
            ss << "NOTE: Provided MPI thread safety level ("
                           << mtn(thread_provided) << ") is stricter than requested "
                           << mtn(thread_required) << "). This might reduce multi-node "
                           << "communication performance.";
            WarpX::GetInstance().RecordWarning("MPI", ss.str());
        }
#endif
    }

} // namespace utils
