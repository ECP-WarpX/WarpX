#include <SpectralFieldDataHankel.H>

#include <WarpX.H>

using namespace amrex;

/* \brief Initialize fields in spectral space, and FFT plans */
SpectralFieldDataHankel::SpectralFieldDataHankel (const BoxArray& realspace_ba,
                                                  const SpectralHankelKSpace& k_space,
                                                  const DistributionMapping& dm,
                                                  const int n_field_required,
                                                  const int n_rz_azimuthal_modes)
    : n_rz_azimuthal_modes(n_rz_azimuthal_modes)
{
    const BoxArray& spectralspace_ba = k_space.spectralspace_ba;

    // Allocate the arrays that contain the fields in spectral space
    // (one component per field)
    fields = SpectralField(spectralspace_ba, dm, n_rz_azimuthal_modes*n_field_required, 0);

    // Allocate temporary arrays - in real space and spectral space
    // These arrays will store the data just before/after the z FFT
    tempHTransformed = SpectralField(realspace_ba, dm, n_rz_azimuthal_modes, 0);
    tmpSpectralField = SpectralField(spectralspace_ba, dm, n_rz_azimuthal_modes, 0);

    // By default, we assume the z FFT is done from/to a nodal grid in real space
    // It the FFT is performed from/to a cell-centered grid in real space,
    // a correcting "shift" factor must be applied in spectral space.
    zshift_FFTfromCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformFromCellCentered);
    zshift_FFTtoCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformToCellCentered);

    // Allocate and initialize the FFT plans and Hankel transformer
    forward_plan = FFTplans(spectralspace_ba, dm);
    backward_plan = FFTplans(spectralspace_ba, dm);
    hankeltransformer = HankelTransformer(spectralspace_ba, dm);

    // Loop over boxes and allocate the corresponding plan
    // for each box owned by the local MPI proc
    for (MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi){
        IntVect grid_size = realspace_ba[mfi].length();
#ifdef AMREX_USE_GPU
        NOT IMPLEMENTED
        // Create cuFFT plans
        // Creating 3D plan for real to complex -- double precision
        // Assuming CUDA is used for programming GPU
        // Note that D2Z is inherently forward plan
        // and  Z2D is inherently backward plan
        cufftResult result;
#if (AMREX_SPACEDIM == 3)
        result = cufftPlan3d(&forward_plan[mfi], grid_size[2],
                             grid_size[1],grid_size[0], CUFFT_D2Z);
        if (result != CUFFT_SUCCESS) {
           amrex::Print() << " cufftplan3d forward failed! \n";
        }

        result = cufftPlan3d(&backward_plan[mfi], grid_size[2],
                             grid_size[1], grid_size[0], CUFFT_Z2D);
        if (result != CUFFT_SUCCESS) {
           amrex::Print() << " cufftplan3d backward failed! \n";
        }
#else
        result = cufftPlan2d(&forward_plan[mfi], grid_size[1],
                              grid_size[0], CUFFT_D2Z);
        if (result != CUFFT_SUCCESS) {
           amrex::Print() << " cufftplan2d forward failed! \n";
        }

        result = cufftPlan2d(&backward_plan[mfi], grid_size[1],
                               grid_size[0], CUFFT_Z2D);
        if (result != CUFFT_SUCCESS) {
           amrex::Print() << " cufftplan2d backward failed! \n";
        }
#endif

#else
        // Create FFTW plans
        fftw_iodim dims[1];
        fftw_iodim howmany_dims[2];
        dims[0].n = grid_size[1];
        dims[0].is = grid_size[0];
        dims[0].os = grid_size[0];
        howmany_dims[0].n = n_rz_azimuthal_modes;
        howmany_dims[0].is = grid_size[0]*grid_size[1];
        howmany_dims[0].os = grid_size[0]*grid_size[1];
        howmany_dims[1].n = grid_size[0];
        howmany_dims[1].is = 1;
        howmany_dims[1].os = 1;
        forward_plan[mfi] =
            // Note that AMReX FAB are Fortran-order
            fftw_plan_guru_dft(1, // int rank
                               dims,
                               2, // int howmany_rank,
                               howmany_dims,
                               reinterpret_cast<fftw_complex*>(tempHTransformed[mfi].dataPtr()), // fftw_complex *in
                               reinterpret_cast<fftw_complex*>(tmpSpectralField[mfi].dataPtr()), // fftw_complex *out
                               FFTW_FORWARD, // int sign
                               FFTW_ESTIMATE); // unsigned flags
        backward_plan[mfi] =
            fftw_plan_guru_dft(1, // int rank
                               dims,
                               2, // int howmany_rank,
                               howmany_dims,
                               reinterpret_cast<fftw_complex*>(tmpSpectralField[mfi].dataPtr()), // fftw_complex *in
                               reinterpret_cast<fftw_complex*>(tempHTransformed[mfi].dataPtr()), // fftw_complex *out
                               FFTW_BACKWARD, // int sign
                               FFTW_ESTIMATE); // unsigned flags

        std::array<Real,3> xmax = WarpX::UpperCorner(mfi.tilebox(), 0); // lev=0 is explicit
        hankeltransformer[mfi] = SpectralHankelTransformer(grid_size[0], n_rz_azimuthal_modes, xmax[2]);
#endif
    }
}


SpectralFieldDataHankel::~SpectralFieldDataHankel()
{
    if (fields.size() > 0){
        for (MFIter mfi(fields); mfi.isValid(); ++mfi){
#ifdef AMREX_USE_GPU
            // Destroy cuFFT plans
            cufftDestroy(forward_plan[mfi]);
            cufftDestroy(backward_plan[mfi]);
#else
            // Destroy FFTW plans
            fftw_destroy_plan(forward_plan[mfi]);
            fftw_destroy_plan(backward_plan[mfi]);
#endif
        }
    }
}

/* \brief Z Transform the FAB to spectral space,
 *  and store the corresponding result internally
 *  (in the spectral field specified by `field_index`) */
void
SpectralFieldDataHankel::FABZForwardTransform (MFIter const & mfi,
                                               amrex::FArrayBox const & tempHTransformedSplit,
                                               const int field_index, const bool is_nodal_z)
{
    const int ncomp = 2*n_rz_azimuthal_modes - 1;

    // Copy the split complex to the interleaved complex
    const Box realspace_bx = tempHTransformed[mfi].box();
    Array4<const Real> split_arr = tempHTransformedSplit.array();
    Array4<Complex> complex_arr = tempHTransformed[mfi].array();
    // mode 0
    {
        ParallelFor(realspace_bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            complex_arr(i,j,k,0) = {split_arr(i,j,k,0), 0.};
        });
    }
    // Can this outer for loop be put inside the ParallelFor?
    for (int mode=1 ; mode < n_rz_azimuthal_modes ; mode++) {
        int icomp = 2*mode - 1;
        ParallelFor(realspace_bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            complex_arr(i,j,k,mode) = {split_arr(i,j,k,icomp), split_arr(i,j,k,icomp+1)};
        });
    }

    // Perform Fourier transform from `tmpRealField` to `tmpSpectralField`
#ifdef AMREX_USE_GPU
    NOT IMPLEMENTED
    // Perform Fast Fourier Transform on GPU using cuFFT
    // make sure that this is done on the same
    // GPU stream as the above copy
    cufftResult result;
    cudaStream_t stream = amrex::Gpu::Device::cudaStream();
    cufftSetStream(forward_plan[mfi], stream);
    result = cufftExecD2Z(forward_plan[mfi],
                          tmpRealField[mfi].dataPtr(),
                          reinterpret_cast<cuDoubleComplex*>(
                          tmpSpectralField[mfi].dataPtr()));
    if (result != CUFFT_SUCCESS) {
       amrex::Print() << " forward transform using cufftExecD2Z failed ! \n";
    }
#else
    fftw_execute(forward_plan[mfi]);
#endif

    // Copy the spectral-space field `tmpSpectralField` to the appropriate
    // index of the FabArray `fields` (specified by `field_index`)
    // and apply correcting shift factor if the real space data comes
    // from a cell-centered grid in real space instead of a nodal grid.
    Array4<const Complex> tmp_arr = tmpSpectralField[mfi].array();
    Array4<Complex> fields_arr = SpectralFieldDataHankel::fields[mfi].array();
    const Complex* zshift_arr = zshift_FFTfromCell[mfi].dataPtr();

    // Loop over indices within one box, all components
    const Box spectralspace_bx = tmpSpectralField[mfi].box();
    const int findex = field_index*n_rz_azimuthal_modes;
    ParallelFor(spectralspace_bx, ncomp,
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        Complex spectral_field_value = tmp_arr(i,j,k,n);
        // Apply proper shift
        if (is_nodal_z==false) spectral_field_value *= zshift_arr[j];
        // Copy field into the right index
        fields_arr(i,j,k,findex+n) = spectral_field_value;
    });
}

/* \brief Backward Z Transform the data from the fields
 * (in the spectral field specified by `field_index`)
 * to physical space, and return the resulting FArrayBox */
amrex::FArrayBox
SpectralFieldDataHankel::FABZBackwardTransform (MFIter const & mfi,
                                                const int field_index, const bool is_nodal_z)
{
    const int ncomp = 2*n_rz_azimuthal_modes - 1;

    // Copy the spectral-space field from the appropriate index of the FabArray
    // `fields` (specified by `field_index`) to field `tmpSpectralField`
    // and apply correcting shift factor if the real space data is on
    // a cell-centered grid in real space instead of a nodal grid.
    Array4<const Complex> fields_arr = SpectralFieldDataHankel::fields[mfi].array();
    Array4<Complex> tmp_arr = tmpSpectralField[mfi].array();
    const Complex* zshift_arr = zshift_FFTtoCell[mfi].dataPtr();

    // Loop over indices within one box, all components
    const Box spectralspace_bx = tmpSpectralField[mfi].box();
    const int findex = field_index*n_rz_azimuthal_modes;
    ParallelFor(spectralspace_bx, ncomp,
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        Complex spectral_field_value = fields_arr(i,j,k,findex+n);
        // Apply proper shift
        if (is_nodal_z==false) spectral_field_value *= zshift_arr[j];
        // Copy field into the right index
        tmp_arr(i,j,k,n) = spectral_field_value;
    });

    // Perform Fourier transform from `tmpSpectralField` to `tmpRealField`
#ifdef AMREX_USE_GPU
    NOT IMPLEMENTED
    // Perform Fast Fourier Transform on GPU using cuFFT.
    // make sure that this is done on the same
    // GPU stream as the above copy
    cufftResult result;
    cudaStream_t stream = amrex::Gpu::Device::cudaStream();
    cufftSetStream(backward_plan[mfi], stream);
    result = cufftExecZ2D(backward_plan[mfi],
                          reinterpret_cast<cuDoubleComplex*>(
                          tmpSpectralField[mfi].dataPtr()),
                          tmpRealField[mfi].dataPtr());
    if (result != CUFFT_SUCCESS) {
       amrex::Print() << " Backward transform using cufftexecZ2D failed! \n";
    }
#else
    fftw_execute(backward_plan[mfi]);
#endif

    // Copy the interleaved complex to the split complex
    Box realspace_bx = tempHTransformed[mfi].box();
    realspace_bx.enclosedCells(); // Discard last point in nodal direction
    amrex::FArrayBox tempHTransformedSplit(realspace_bx, ncomp);
    AMREX_ALWAYS_ASSERT(realspace_bx == tempHTransformed[mfi].box());
    Array4<const Complex> complex_arr = tempHTransformed[mfi].array();
    Array4<Real> split_arr = tempHTransformedSplit.array();
    // mode 0
    {
        ParallelFor(realspace_bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            split_arr(i,j,k,0) = complex_arr(i,j,k,0).real();
        });
    }
    // Can this outer for loop be put inside the ParallelFor?
    for (int mode=1 ; mode < n_rz_azimuthal_modes ; mode++) {
        int icomp = 2*mode - 1;
        ParallelFor(realspace_bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            split_arr(i,j,k,icomp) = complex_arr(i,j,k,mode).real();
            split_arr(i,j,k,icomp+1) = complex_arr(i,j,k,mode).imag();
        });
    }

    return tempHTransformedSplit;

}

/* \brief Transform the component `i_comp` of MultiFab `mf`
 *  to spectral space, and store the corresponding result internally
 *  (in the spectral field specified by `field_index`) */
void
SpectralFieldDataHankel::ForwardTransform (const MultiFab& mf, const int field_index,
                                           const int i_comp)
{
    // Check field index type, in order to apply proper shift in spectral space
    // Only cell centered in r is supported.
    const bool is_nodal_z = mf.is_nodal(1);

    const int ncomp = 2*n_rz_azimuthal_modes - 1;

    // Loop over boxes
    for (MFIter mfi(mf); mfi.isValid(); ++mfi){

        // Perform the Hankel transform first
        Box realspace_bx = mf[mfi].box(); // Copy the box
        realspace_bx.enclosedCells(); // Discard last point in nodal direction
        amrex::FArrayBox tempHTransformedSplit(realspace_bx, ncomp);
        amrex::FArrayBox mf_comp(mf[mfi], amrex::make_alias, i_comp*ncomp, ncomp);
        hankeltransformer[mfi].PhysicalToSpectral_Scalar(realspace_bx, mf_comp, tempHTransformedSplit);

        FABZForwardTransform(mfi, tempHTransformedSplit, field_index, is_nodal_z);

    }
}

/* \brief Transform the coupled components of MultiFabs `mf_r` and `mf_t`
 *  to spectral space, and store the corresponding result internally
 *  (in the spectral fields specified by `field_index_r` and `field_index_r`) */
void
SpectralFieldDataHankel::ForwardTransform (const MultiFab& mf_r, const int field_index_r,
                                           const MultiFab& mf_t, const int field_index_t)
{
    // Check field index type, in order to apply proper shift in spectral space
    // Only cell centered in r is supported.
    const bool is_nodal_z = mf_r.is_nodal(1);

    const int ncomp = 2*n_rz_azimuthal_modes - 1;

    // Loop over boxes
    for (MFIter mfi(mf_r); mfi.isValid(); ++mfi){

        // Perform the Hankel transform first
        Box realspace_bx = mf_r[mfi].box(); // Copy the box
        realspace_bx.enclosedCells(); // Discard last point in nodal direction
        amrex::FArrayBox tempHTransformedSplit_p(realspace_bx, ncomp);
        amrex::FArrayBox tempHTransformedSplit_m(realspace_bx, ncomp);
        hankeltransformer[mfi].PhysicalToSpectral_Vector(realspace_bx,
                                                         mf_r[mfi], mf_t[mfi],
                                                         tempHTransformedSplit_p, tempHTransformedSplit_m);

        FABZForwardTransform(mfi, tempHTransformedSplit_p, field_index_r, is_nodal_z);
        FABZForwardTransform(mfi, tempHTransformedSplit_m, field_index_t, is_nodal_z);

    }
}

/* \brief Transform spectral field specified by `field_index` back to
 * real space, and store it in the component `i_comp` of `mf` */
void
SpectralFieldDataHankel::BackwardTransform (MultiFab& mf, const int field_index,
                                            const int i_comp)
{
    // Check field index type, in order to apply proper shift in spectral space
    const bool is_nodal_z = mf.is_nodal(1);

    const int ncomp = 2*n_rz_azimuthal_modes - 1;

    // Loop over boxes
    for (MFIter mfi(mf); mfi.isValid(); ++mfi){

        amrex::FArrayBox tempHTransformedSplit = FABZBackwardTransform(mfi, field_index, is_nodal_z);

        // Perform the Hankel inverse transform last
        const Box realspace_bx = mf[mfi].box(); // Copy the box
        amrex::FArrayBox mf_comp(mf[mfi], amrex::make_alias, i_comp*ncomp, ncomp);
        hankeltransformer[mfi].SpectralToPhysical_Scalar(realspace_bx, tempHTransformedSplit, mf_comp);

        // Copy the temporary field `tmpRealField` to the real-space field `mf`
        // (only in the valid cells ; not in the guard cells)
        // Normalize (divide by 1/N) since the FFT+IFFT results in a factor N
        /* ?????
        {
            Array4<Real> mf_arr = mf[mfi].array();
            Array4<const Real> tmp_arr = tmpRealField[mfi].array();
            // Normalization: divide by the number of points in realspace
            // (includes the guard cells)
            const Box realspace_bx = tmpRealField[mfi].box();
            const Real inv_N = 1./realspace_bx.numPts();

            ParallelFor(mfi.validbox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Copy and normalize field
                mf_arr(i,j,k,i_comp) = inv_N*tmp_arr(i,j,k);
            });
        }
        */
    }
}

/* \brief Transform spectral fields specified by `field_index_r` and
 * `field_index_t` back to real space, and store them in `mf_r` and `mf_t` */
void
SpectralFieldDataHankel::BackwardTransform (MultiFab& mf_r, const int field_index_r,
                                            MultiFab& mf_t, const int field_index_t)
{
    // Check field index type, in order to apply proper shift in spectral space
    const bool is_nodal_z = mf_r.is_nodal(1);

    const int ncomp = 2*n_rz_azimuthal_modes - 1;

    // Loop over boxes
    for (MFIter mfi(mf_r); mfi.isValid(); ++mfi){

        amrex::FArrayBox tempHTransformedSplit_p = FABZBackwardTransform(mfi, field_index_r, is_nodal_z);
        amrex::FArrayBox tempHTransformedSplit_m = FABZBackwardTransform(mfi, field_index_t, is_nodal_z);

        // Perform the Hankel inverse transform last
        const Box realspace_bx = mf_r[mfi].box(); // Copy the box
        hankeltransformer[mfi].SpectralToPhysical_Vector(realspace_bx,
                                                         tempHTransformedSplit_p, tempHTransformedSplit_m,
                                                         mf_r[mfi], mf_t[mfi]);

        // Copy the temporary field `tmpRealField` to the real-space field `mf`
        // (only in the valid cells ; not in the guard cells)
        // Normalize (divide by 1/N) since the FFT+IFFT results in a factor N
        /* ?????
        {
            Array4<Real> mf_arr = mf[mfi].array();
            Array4<const Real> tmp_arr = tmpRealField[mfi].array();
            // Normalization: divide by the number of points in realspace
            // (includes the guard cells)
            const Box realspace_bx = tmpRealField[mfi].box();
            const Real inv_N = 1./realspace_bx.numPts();

            ParallelFor(mfi.validbox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Copy and normalize field
                mf_arr(i,j,k,i_comp) = inv_N*tmp_arr(i,j,k);
            });
        }
        */
    }
}
