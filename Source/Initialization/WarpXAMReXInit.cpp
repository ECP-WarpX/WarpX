/* Copyright 2020 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Initialization/WarpXAMReXInit.H"

#include <AMReX_ParmParse.H>


namespace {
    /** Overwrite defaults in AMReX Inputs
     *
     * This overwrites defaults in amrex::ParamParse for inputs.
     */
    void
    overwrite_amrex_parser_defaults()
    {
        amrex::ParmParse pp("amrex");

        // https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters
        bool abort_on_out_of_gpu_memory = true; // AMReX' default: false
        pp.query("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);
        pp.add("abort_on_out_of_gpu_memory", abort_on_out_of_gpu_memory);
    }
}

amrex::AMReX*
warpx_amrex_init(int& argc, char**& argv)
{
    // note: AMReX defines a placeholder/"mock-up" for MPI_COMM_WORLD in serial builds
    bool const build_parm_parse = true; // AMReX' default
    MPI_Comm const mpi_comm = MPI_COMM_WORLD; // AMReX' default

    return amrex::Initialize(
        argc,
        argv,
        build_parm_parse,
        mpi_comm,
        overwrite_amrex_parser_defaults
    );
}
