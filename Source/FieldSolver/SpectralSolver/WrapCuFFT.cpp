#include "AnyFFT.H"

namespace AnyFFT
{
#ifdef AMREX_USE_FLOAT
    const auto mapdir std::map<direction, cufftType>{
        {R2C, CUFFT_R2C},
        {C2R, CUFFT_C2R}};
#else
    const auto mapdir std::map<direction, cufftType>{
        {R2C, CUFFT_D2Z},
        {C2R, CUFFT_Z2D}};
#endif

    std::string cufftErrorToString (const cufftResult& err);

    FFTplan CreatePlan(int nx, int ny, int nz,
                       amrex::Real* real_array, Complex* complex_array, direction dir)
    {
        FFTplan fft_plan;

#if (AMREX_SPACEDIM == 3)
        cufftResult result = cufftPlan3d( &(fft_plan.m_plan), nz, ny, nx, mapdir[direction]);
#else
        cufftResult result = cufftPlan2d( &(fft_plan.m_plan), nz, ny, nx, mapdir[direction]);
#endif
        if ( result != CUFFT_SUCCESS ) {
            amrex::Print() << " cufftplan failed! Error: " <<
                cufftErrorToString(result) << "\n";
        }

        fft_plan.m_real_array = real_array;
        fft_plan.m_complex_array = complex_array;
        fft_plan.m_dir = dir;

        return fft_plan;
    }

    void DestroyPlan(FFTplan& fft_plan)
    {
        cufftDestroy( fft_plan.m_plan );
    }

    void Execute(FFTplan& fft_plan){
        // make sure that this is done on the same GPU stream as the above copy
        cudaStream_t stream = amrex::Gpu::Device::cudaStream();
        cufftSetStream ( fft_plan.m_plan, stream);
        cufftResult result;
        if (fft_plan.m_dir == direction::R2C){
#  ifdef AMREX_USE_FLOAT
            result = cufftExecR2C(fft_plan.m_plan, fft_plan.m_real_array, fft_plan.m_complex_array);
#  else
            result = cufftExecD2Z(fft_plan.m_plan, fft_plan.m_real_array, fft_plan.m_complex_array);
#  endif
        } else if (fft_plan.m_dir == direction::C2R){
#  ifdef AMREX_USE_FLOAT
            result = cufftExecC2R(fft_plan.m_plan, fft_plan.m_complex_array, fft_plan.m_real_array);
#  else
            result = cufftExecZ2D(fft_plan.m_plan, fft_plan.m_complex_array, fft_plan.m_real_array);
#  endif
        } else {
            amrex::Abort("direction must be AnyFFT::direction::R2C or AnyFFT::direction::C2R");
        }
        if ( result != CUFFT_SUCCESS ) {
            amrex::Print() << " forward transform using cufftExec failed ! Error: " <<
                cufftErrorToString(result) << "\n";
        }
    }
    
    /** \brief This method converts a cufftResult
     * into the corresponding string
     *
     * @param[in] err a cufftResult
     * @return an std::string
     */
    std::string cufftErrorToString (const cufftResult& err)
    {
        const auto res2string = std::map<cufftResult, std::string>{
            {CUFFT_SUCCESS, "CUFFT_SUCCESS"},
            {CUFFT_INVALID_PLAN,"CUFFT_INVALID_PLAN"},
            {CUFFT_ALLOC_FAILED,"CUFFT_ALLOC_FAILED"},
            {CUFFT_INVALID_TYPE,"CUFFT_INVALID_TYPE"},
            {CUFFT_INVALID_VALUE,"CUFFT_INVALID_VALUE"},
            {CUFFT_INTERNAL_ERROR,"CUFFT_INTERNAL_ERROR"},
            {CUFFT_EXEC_FAILED,"CUFFT_EXEC_FAILED"},
            {CUFFT_SETUP_FAILED,"CUFFT_SETUP_FAILED"},
            {CUFFT_INVALID_SIZE,"CUFFT_INVALID_SIZE"},
            {CUFFT_UNALIGNED_DATA,"CUFFT_UNALIGNED_DATA"}};

        const auto it = res2string.find(err);
        if(it != res2string.end()){
            return it->second;
        }
        else{
            return std::to_string(err) +
                " (unknown error code)";
        }
    }
}
