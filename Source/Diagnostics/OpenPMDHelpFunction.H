/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Juliette Pech, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_OPENPMDHELPFUNCTION_H_
#define WARPX_OPENPMDHELPFUNCTION_H_

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif
#include <string>


std::string
WarpXOpenPMDFileType ();

#endif  // WARPX_OPENPMDHELPFUNCTION_H_
