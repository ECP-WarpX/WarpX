#ifndef WRAPFFTW_H_
#define WRAPFFTW_H_

#include <fftw3.h>

namespace AnyFFT
{
#ifdef AMREX_USE_FLOAT
    using Complex = fftwf_complex; /**< Complex type for FFTW in single precision */
    /** Define the FFTplans type, which holds one fft plan per box
     * (plans are only initialized for the boxes that are owned by the local MPI rank). */
    using VendorFFTPlan = fftwf_plan;
#else
    /** Define the FFTplans type, which holds one fft plan per box
     * (plans are only initialized for the boxes that are owned by the local MPI rank). */
    using Complex = fftw_complex; /**< Complex type for FFTW in double precision */
    using VendorFFTPlan = fftw_plan;
#endif

    /**
     * This struct contains the vendor FFT plan and additional metadata
     */
    struct FFTplan
    {
        amrex::Real* m_real_array; /**< pointer to real array */
        Complex* m_complex_array; /**< pointer to complex array */
        VendorFFTPlan m_plan; /**< FFTW plan */
        direction m_dir;  /**< direction (C2R or R2C) */
    };

    /** Collection of FFT plans, one FFTplan per box */
    using FFTplans = amrex::LayoutData<FFTplan>;

    /** \brief create FFTW plan
     * \param[in] nx array size in 1st direction
     * \param[in] ny array size in 2nd direction
     * \param[in] nz array size in 3rd direction
     * \param[out] real_array Real array from/to where R2C/C2R FFT is performed
     * \param[out] complex_array Complex array to/from where R2C/C2R FFT is performed
     * \param[in] dir direction, either R2C or C2R
     */
    FFTplan CreatePlan(int nx, int ny, int nz,
                       amrex::Real* real_array, Complex* complex_array, direction dir);

    /** \brief Destroy FFTW plan
     * \param[out] fft_plan plan to destroy
     */
    void DestroyPlan(FFTplan& fft_plan);

    /** \brief Perform FFT with FFTW
     * \param[out] fft_plan plan for which the FFT is performed
     */
    void Execute(FFTplan& fft_plan);
}

#endif // WRAPFFTW_H_
