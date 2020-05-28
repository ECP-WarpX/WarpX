#ifndef WRAPCUFFT_H_
#define WRAPCUFFT_H_

#include <cufft.h>

namespace AnyFFT
{
    /** Define the FFTplans type, which holds one fft plan per box
     * (plans are only initialized for the boxes that are owned by the local MPI rank). */
    using VendorFFTPlan = cufftHandle;
#ifdef AMREX_USE_FLOAT
    using Complex = cuComplex; /**< Complex type for cuFFT in single precision */
#else
    using Complex = cuDoubleComplex; /**< Complex type for cuFFT in double precision */
#endif

    /**
     * This struct contains the vendor FFT plan and additional metadata
     */
    struct FFTplan
    {
        amrex::Real* m_real_array; /**< pointer to real array */
        Complex* m_complex_array; /**< pointer to complex array */
        VendorFFTPlan m_plan; /**< cuFFT plan */
        direction m_dir;  /**< direction (C2R or R2C) */
    };

    /** Collection of FFT plans, one FFTplan per box */
    using FFTplans = amrex::LayoutData<FFTplan>;

    /** \brief create cuFFT plan assuming CUDA is used for programming GPU
     * \param[in] real_size size of the real array
     * \param[out] real_array Real array from/to where R2C/C2R FFT is performed
     * \param[out] complex_array Complex array to/from where R2C/C2R FFT is performed
     * \param[in] dir direction, either R2C or C2R
     */
    FFTplan CreatePlan(IntVect real_size, amrex::Real* real_array,
                       Complex* complex_array, direction dir);

    /** \brief Destroy cuFFT plan
     * \param[out] fft_plan plan to destroy
     */
    void DestroyPlan(FFTplan& fft_plan);

    /** \brief Perform FFT with cuFFT
     * \param[out] fft_plan plan for which the FFT is performed
     */
    void Execute(FFTplan& fft_plan);
}

#endif
