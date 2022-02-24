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
    /** Return the required MPI threading
     *
     * @return the MPI_THREAD_* level required for MPI_Init_thread
     */
    int
    warpx_mpi_thread_required ();

    /** Initialize MPI
     *
     * @return pair(required, provided) of MPI thread level from MPI_Init_thread
     */
    std::pair< int, int >
    warpx_mpi_init (int argc, char* argv[]);

    /** Check if the requested MPI thread level is valid
     *
     * Prints warnings and notes otherwise.
     */
    void
    warpx_check_mpi_thread_level ();

} // namespace utils

#endif // WARPX_MPI_INIT_HELPERS_H_
