/* Copyright 2020 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_AMREX_INIT_H_
#define WARPX_AMREX_INIT_H_

#include <AMReX.H>

/** Call amrex::Initialize
 *
 * This function calls amrex::Initialize and overwrites AMReX' defaults.
 *
 * @param[in] argc number of arguments from main()
 * @param[in] argv argument strings from main()
 * @returns pointer to an AMReX* object, forwarded from amrex::Initialize
 */
amrex::AMReX*
warpx_amrex_init(int& argc, char**& argv);

#endif
