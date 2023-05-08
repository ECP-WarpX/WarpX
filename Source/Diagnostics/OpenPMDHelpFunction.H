//
// Created by juliette on 03/05/23.
//
#ifndef WARPX_OPENPMDHELPFUNCTION_H_
#define WARPX_OPENPMDHELPFUNCTION_H_

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif
#include <string>


std::string
WarpXOpenPMDFileType ();

#endif  // WARPX_OPENPMDHELPFUNCTION_H_
