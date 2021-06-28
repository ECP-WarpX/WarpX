/* Copyright 2020 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MPI_INIT_HELPERS_H_
#define WARPX_MPI_INIT_HELPERS_H_

#include <utility>

namespace utils
{
    /** Initialize MPI
     *
     * @return pair(required, provided) of MPI thread level from MPI_Init_thread
     */
    std::pair< int, int >
    warpx_mpi_init (int argc, char* argv[]);

    /** Check if the requested MPI thread level is valid
     *
     * Prints warnings and notes otherwise.
     *
     * @param mpi_thread_levels pair(required, provided) of MPI thread level from MPI_Init_thread
     */
    void
    warpx_check_mpi_thread_level (std::pair< int, int > const mpi_thread_levels);

} // namespace utils

#endif // WARPX_MPI_INIT_HELPERS_H_
