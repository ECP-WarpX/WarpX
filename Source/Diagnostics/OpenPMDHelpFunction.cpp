/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Juliette Pech, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "OpenPMDHelpFunction.H"
#include "Utils/TextMsg.H"

std::string
WarpXOpenPMDFileType ()
{
    std::string openPMDFileType;
#ifdef WARPX_USE_OPENPMD
#if openPMD_HAVE_ADIOS2==1
    openPMDFileType = "bp";
#elif openPMD_HAVE_ADIOS1==1
    openPMDFileType = "bp";
#elif openPMD_HAVE_HDF5==1
    openPMDFileType = "h5";
#else
    openPMDFileType = "json";
#endif
#else
    WARPX_ABORT_WITH_MESSAGE("openPMD-api cannot be used!");
#endif // WARPX_USE_OPENPMD
    return openPMDFileType;
}
