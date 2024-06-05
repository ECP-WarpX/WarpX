/* Copyright 2020 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Initialization/WarpXAMReXInit.H"

#include "Utils/TextMsg.H"

#include <AMReX.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TinyProfiler.H>

#include <string>

namespace {
    /** Overwrite defaults in AMReX Inputs
     *
     * This overwrites defaults in amrex::ParmParse for inputs.
     */
    void
    overwrite_amrex_parser_defaults ()
    {
        amrex::ParmParse pp_amrex("amrex");

        // https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters
        bool abort_on_out_of_gpu_memory = true; // AMReX' default: false
        pp_amrex.queryAdd("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);

        bool the_arena_is_managed = false; // AMReX' default: true
        pp_amrex.queryAdd("the_arena_is_managed", the_arena_is_managed);

        // https://amrex-codes.github.io/amrex/docs_html/InputsComputeBackends.html
        std::string omp_threads = "nosmt"; // AMReX' default: system
        pp_amrex.queryAdd("omp_threads", omp_threads);

        // Work-around:
        // If warpx.numprocs is used for the domain decomposition, we will not use blocking factor
        // to generate grids. Nonetheless, AMReX has asserts in place that validate that the
        // number of cells is a multiple of blocking factor. We set the blocking factor to 1 so those
        // AMReX asserts will always pass.
        const amrex::ParmParse pp_warpx("warpx");
        if (pp_warpx.contains("numprocs"))
        {
            amrex::ParmParse pp_amr("amr");
            pp_amr.add("blocking_factor", 1);
        }

        //See https://github.com/AMReX-Codes/amrex/pull/3763
#ifdef AMREX_USE_GPU
        bool warpx_do_device_synchronize = true;
#else
        bool warpx_do_device_synchronize = false;
#endif
        pp_warpx.query("do_device_synchronize", warpx_do_device_synchronize);
        bool do_device_synchronize = warpx_do_device_synchronize;
        amrex::ParmParse pp_tiny_profiler("tiny_profiler");
        if (pp_tiny_profiler.queryAdd("device_synchronize_around_region", do_device_synchronize) )
        {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                do_device_synchronize == warpx_do_device_synchronize,
                "tiny_profiler.device_synchronize_around_region overrides warpx.do_device_synchronize.");
        }

        // Here we override the default tiling option for particles, which is always
        // "false" in AMReX, to "false" if compiling for GPU execution and "true"
        // if compiling for CPU.
        {
            amrex::ParmParse pp_particles("particles");
#ifdef AMREX_USE_GPU
            bool do_tiling = false; // By default, tiling is off on GPU
#else
            bool do_tiling = true;
#endif
            pp_particles.queryAdd("do_tiling", do_tiling);
        }
    }
}

namespace warpx::initialization
{

    amrex::AMReX*
    amrex_init (int& argc, char**& argv, bool build_parm_parse)
    {
        return amrex::Initialize(
            argc,
            argv,
            build_parm_parse,
            MPI_COMM_WORLD,
            ::overwrite_amrex_parser_defaults
        );
    }

}
