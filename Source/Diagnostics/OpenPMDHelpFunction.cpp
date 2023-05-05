//
// Created by juliette on 03/05/23.
//
#include "OpenPMDHelpFunction.H"
#include "WarpXOpenPMD.H"

#ifdef WARPX_USE_OPENPMD
std::string
WarpXOpenPMDFileType ()
{
    std::string openPMDFileType;
#if openPMD_HAVE_ADIOS2==1
    openPMDFileType = "bp";
#elif openPMD_HAVE_ADIOS1==1
    openPMDFileType = "bp";
#elif openPMD_HAVE_HDF5==1
    openPMDFileType = "h5";
#else
    openPMDFileType = "json";
#endif
    return openPMDFileType;
}
#endif // WARPX_USE_OPENPMD
