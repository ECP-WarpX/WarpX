/* Copyright 2020 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Initialization/WarpXAMReXInit.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <memory>

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

        // Work-around:
        // If warpx.numprocs is used for the domain decomposition, we will not use blocking factor
        // to generate grids. Nonetheless, AMReX has asserts in place that validate that the
        // number of cells is a multiple of blocking factor. We set the blocking factor to 1 so those
        // AMReX asserts will always pass.
        amrex::ParmParse pp_warpx("warpx");
        if (pp_warpx.contains("numprocs"))
        {
            amrex::ParmParse pp_amr("amr");
            pp_amr.add("blocking_factor", 1);
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

amrex::AMReX*
warpx_amrex_init (int& argc, char**& argv, bool const build_parm_parse, MPI_Comm const mpi_comm)
{
    return amrex::Initialize(
        argc,
        argv,
        build_parm_parse,
        mpi_comm,
        overwrite_amrex_parser_defaults
    );
}
