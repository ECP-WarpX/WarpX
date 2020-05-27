#include "AnyFFT.H"

namespace AnyFFT
{
FFTplan CreatePlan(int nx, int ny, int nz,
        amrex::Real* real_array, PrecisionComplex* complex_array, direction dir)
    {
        FFTplan fft_plan;
        
        if (dir == direction::R2C){
        // Create FFTW plans
        fft_plan.m_plan =
            // Swap dimensions: AMReX FAB are Fortran-order but FFTW is C-order
#if (AMREX_SPACEDIM == 3)
#  ifdef AMREX_USE_FLOAT
            fftwf_plan_dft_r2c_3d( nz, ny, nx, real_array, complex_array, FFTW_ESTIMATE);
#  else
            fftw_plan_dft_r2c_3d( nz, ny, nx, real_array, complex_array, FFTW_ESTIMATE);
#  endif
#else
#  ifdef AMREX_USE_FLOAT
            fftwf_plan_dft_r2c_2d( ny, nx, real_array, complex_array, FFTW_ESTIMATE);
#  else
            fftw_plan_dft_r2c_2d( ny, nx, real_array, complex_array, FFTW_ESTIMATE);
#  endif
#endif
        } else if (dir == direction::C2R){
            // Swap dimensions: AMReX FAB are Fortran-order but FFTW is C-order
#if (AMREX_SPACEDIM == 3)
#  ifdef AMREX_USE_FLOAT
            fftwf_plan_dft_c2r_3d( nz, ny, nx, complex_array, real_array, FFTW_ESTIMATE);
#  else
            fftw_plan_dft_c2r_3d( nz, ny, nx, complex_array, real_array, FFTW_ESTIMATE);
#  endif
#else
#  ifdef AMREX_USE_FLOAT
            fftwf_plan_dft_c2r_2d( ny, nx, complex_array, real_array, FFTW_ESTIMATE);
#  else
            fftw_plan_dft_c2r_2d( ny, nx, complex_array, real_array, FFTW_ESTIMATE);
#  endif
#endif
        }
        fft_plan.m_real_array = real_array;
        fft_plan.m_complex_array = complex_array;
        fft_plan.m_dir = dir;
        return fft_plan;
    }

    void DestroyPlan(FFTplan fft_plan)
    {
#  ifdef AMREX_USE_FLOAT
        fftwf_destroy_plan( fft_plan.m_plan );
#  else
        fftw_destroy_plan( fft_plan.m_plan );
#  endif        
    }

    void Execute(FFTplan fft_plan){
#  ifdef AMREX_USE_FLOAT
        fftwf_execute( fft_plan.m_plan );
#  else
        fftw_execute( fft_plan.m_plan );
#  endif
    }
}
