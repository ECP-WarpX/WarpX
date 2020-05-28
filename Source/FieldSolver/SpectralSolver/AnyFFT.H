#ifndef ANYFFT_H_
#define ANYFFT_H_

#ifdef AMREX_USE_GPU
#  include <cufft.h>
#else
#  include <fftw3.h>
#endif

#include <AMReX_LayoutData.H>

namespace AnyFFT
{
    // First, define library-independent types (complex, FFT plan)
    
    /** Complex type for FFT, depends on FFT library */
#ifdef AMREX_USE_GPU
#  ifdef AMREX_USE_FLOAT
    using Complex = cuComplex;
#  else
    using Complex = cuDoubleComplex;
#  endif
#else
#  ifdef AMREX_USE_FLOAT
    using Complex = fftwf_complex;
#  else
    using Complex = fftw_complex;
#  endif
#endif
    
    /** Library-dependent FFT plans type, which holds one fft plan per box
     * (plans are only initialized for the boxes that are owned by the local MPI rank). */
#ifdef AMREX_USE_GPU
    using VendorFFTPlan = cufftHandle;
#else
#  ifdef AMREX_USE_FLOAT
    using VendorFFTPlan = fftwf_plan;
#  else
    using VendorFFTPlan = fftw_plan;
#  endif
#endif

    // Second, define library-independent types and functions

    /** Direction in which the FFT is performed. */
    enum struct direction {R2C, C2R};

    /** This struct contains the vendor FFT plan and additional metadata
     */
    struct FFTplan
    {
        amrex::Real* m_real_array; /**< pointer to real array */
        Complex* m_complex_array; /**< pointer to complex array */
        VendorFFTPlan m_plan; /**< Vendor FFT plan */
        direction m_dir;  /**< direction (C2R or R2C) */
    };
    
    /** Collection of FFT plans, one FFTplan per box */
    using FFTplans = amrex::LayoutData<FFTplan>;

    /** \brief create FFT plan for the backend FFT library.
     * \param[in] real_size size of the real array
     * \param[out] real_array Real array from/to where R2C/C2R FFT is performed
     * \param[out] complex_array Complex array to/from where R2C/C2R FFT is performed
     * \param[in] dir direction, either R2C or C2R
     */
    FFTplan CreatePlan(const amrex::IntVect& real_size, amrex::Real * const real_array,
                       Complex * const complex_array, const direction dir);

    /** \brief Destroy library FFT plan.
     * \param[out] fft_plan plan to destroy
     */
    void DestroyPlan(FFTplan& fft_plan);

    /** \brief Perform FFT with backend library.
     * \param[out] fft_plan plan for which the FFT is performed
     */
    void Execute(FFTplan& fft_plan);
}

#endif // ANYFFT_H_
