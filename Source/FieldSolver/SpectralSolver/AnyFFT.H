#ifndef ANYFFT_H_
#define ANYFFT_H_

#include <AMReX_LayoutData.H>

namespace AnyFFT
{
    enum struct direction {R2C, C2R};
}

#ifdef AMREX_USE_GPU
#   include "WrapCuFFT.H"
#else
#   include "WrapFFTW.H"
#endif

#endif
