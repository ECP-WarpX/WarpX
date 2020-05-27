#ifndef WRAPCUFFT_H_
#define WRAPCUFFT_H_

namespace AnyFFT
{    
#ifdef AMREX_USE_FLOAT
    using PrecisionComplex = cuComplex;
#else
    using PrecisionComplex = cuDoubleComplex;
#endif

    using FFTplans = amrex::LayoutData<cufftHandle>;

#ifdef AMREX_USE_FLOAT
    const auto mapdir std::map<direction, cufftType>{
        {R2C, CUFFT_R2C},
        {C2R, CUFFT_C2R}};
#else
    const auto mapdir std::map<direction, cufftType>{
        {R2C, CUFFT_D2Z},
        {C2R, CUFFT_Z2D}};
#endif

    struct FFTPlans {
        amrex::Real* m_real_array;
        PrecisionComplex* m_complex_array;
        amrex::LayoutData<cufftHandle> m_plan;
    }

    CreatePlan(
        FFTPlans& fft_plan, int nx, int ny, int nz,
        amrex::Real* real_array, PrecisionComplex* complex_array, direction dir)
    {
        // Create cuFFT plans
        // Creating 3D plan for real to complex -- double precision
        // Assuming CUDA is used for programming GPU
        // Note that D2Z is inherently forward plan
        // and  Z2D is inherently backward plan
        cufftResult result;
#if (AMREX_SPACEDIM == 3)
        result = cufftPlan3d( &(fft_plan.m_plan), nz, ny, nx, mapdir[direction]);
#else
        result = cufftPlan2d( &(fft_plan.m_plan), nz, ny, nx, mapdir[direction]);
#endif
        if ( result != CUFFT_SUCCESS ) {
            amrex::Print() << " cufftplan failed! Error: " <<
                cufftErrorToString(result) << "\n";
        }
        fft_plan.m_real_array = real_array;
        fft_plan.m_complex_array = complex_array;
    }
}

#endif
