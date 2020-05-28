#include "AnyFFT.H"

namespace AnyFFT
{
#ifdef AMREX_USE_FLOAT
#  if (AMREX_SPACEDIM == 3)
    const auto VendorCreatePlanR2C3D = fftwf_plan_dft_r2c_3d;
    const auto VendorCreatePlanC2R3D = fftwf_plan_dft_c2r_3d;
#  else
    const auto VendorCreatePlanR2C2D = fftwf_plan_dft_r2c_2d;
    const auto VendorCreatePlanC2R2D = fftwf_plan_dft_c2r_2d;
#  endif
#else
#  if (AMREX_SPACEDIM == 3)
    const auto VendorCreatePlanR2C3D = fftw_plan_dft_r2c_3d;
    const auto VendorCreatePlanC2R3D = fftw_plan_dft_c2r_3d;
#  else
    const auto VendorCreatePlanR2C2D = fftw_plan_dft_r2c_2d;
    const auto VendorCreatePlanC2R2D = fftw_plan_dft_c2r_2d;
#  endif
#endif

    FFTplan CreatePlan(const amrex::IntVect& real_size, amrex::Real * const real_array,
                       Complex * const complex_array, const direction dir)
    {
        FFTplan fft_plan;

        if (dir == direction::R2C){
            // Swap dimensions: AMReX FAB are Fortran-order but FFTW is C-order
#if (AMREX_SPACEDIM == 3)
            fft_plan.m_plan = VendorCreatePlanR2C3D(
                real_size[2], real_size[1], real_size[0], real_array, complex_array, FFTW_ESTIMATE);
#else
            fft_plan.m_plan = VendorCreatePlanR2C2D(
                real_size[1], real_size[0], real_array, complex_array, FFTW_ESTIMATE);
#endif
        } else if (dir == direction::C2R){
            // Swap dimensions: AMReX FAB are Fortran-order but FFTW is C-order
#if (AMREX_SPACEDIM == 3)
            fft_plan.m_plan = VendorCreatePlanC2R3D(
                real_size[2], real_size[1], real_size[0], complex_array, real_array, FFTW_ESTIMATE);
#else
            fft_plan.m_plan = VendorCreatePlanC2R2D(
                real_size[1], real_size[0], complex_array, real_array, FFTW_ESTIMATE);
#endif
        }

        fft_plan.m_real_array = real_array;
        fft_plan.m_complex_array = complex_array;
        fft_plan.m_dir = dir;

        return fft_plan;
    }

    void DestroyPlan(FFTplan& fft_plan)
    {
#  ifdef AMREX_USE_FLOAT
        fftwf_destroy_plan( fft_plan.m_plan );
#  else
        fftw_destroy_plan( fft_plan.m_plan );
#  endif
    }

    void Execute(FFTplan& fft_plan){
#  ifdef AMREX_USE_FLOAT
        fftwf_execute( fft_plan.m_plan );
#  else
        fftw_execute( fft_plan.m_plan );
#  endif
    }
}
