#ifndef WRAPFFTW_H_
#define WRAPFFTW_H_

#include <fftw3.h>

namespace AnyFFT
{
#ifdef AMREX_USE_FLOAT
    using PrecisionComplex = fftwf_complex;
    using VendorFFTPlan = fftwf_plan;
#else
    using PrecisionComplex = fftw_complex;
    using VendorFFTPlan = fftw_plan;
#endif

    struct FFTplan
    {
        amrex::Real* m_real_array;
        PrecisionComplex* m_complex_array;
        VendorFFTPlan m_plan;
        direction m_dir;
    };

    using FFTplans = amrex::LayoutData<FFTplan>;

    FFTplan CreatePlan(int nx, int ny, int nz,
                       amrex::Real* real_array, PrecisionComplex* complex_array, direction dir);

    void DestroyPlan(FFTplan fft_plan);

    void Execute(FFTplan fft_plan);
}

#endif
