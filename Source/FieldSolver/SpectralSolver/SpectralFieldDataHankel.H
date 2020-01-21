#ifndef WARPX_SPECTRAL_FIELD_DATA_HANKEL_H_
#define WARPX_SPECTRAL_FIELD_DATA_HANKEL_H_

#include <WarpX_ComplexForFFT.H>
#include <SpectralHankelKSpace.H>
#include <SpectralFieldData.H>
#include <SpectralHankelTransformer.H>
#include <AMReX_MultiFab.H>

/* \brief Class that stores the fields in spectral space, and performs the
 *  Fourier transforms between real space and spectral space
 */
class SpectralFieldDataHankel
{

    // Define the FFTplans type, which holds one fft plan per box
    // (plans are only initialized for the boxes that are owned by
    // the local MPI rank)
#ifdef AMREX_USE_GPU
    using FFTplans = amrex::LayoutData<cufftHandle>;
#else
    using FFTplans = amrex::LayoutData<fftw_plan>;
#endif
    using HankelTransformer = amrex::LayoutData<SpectralHankelTransformer>;

    public:
        SpectralFieldDataHankel(const amrex::BoxArray& realspace_ba,
                                const SpectralHankelKSpace& k_space,
                                const amrex::DistributionMapping& dm,
                                const int n_field_required,
                                const int n_modes);
        SpectralFieldDataHankel() = default; // Default constructor
        SpectralFieldDataHankel& operator=(SpectralFieldDataHankel&& field_data) = default;
        ~SpectralFieldDataHankel();
        void ForwardTransform(const amrex::MultiFab& mf,
                              const int field_index, const int i_comp);
        void ForwardTransform(const amrex::MultiFab& mf_r, const int field_index_r,
                              const amrex::MultiFab& mf_t, const int field_index_t);
        void BackwardTransform(amrex::MultiFab& mf,
                               const int field_index, const int i_comp);
        void BackwardTransform(amrex::MultiFab& mf_r, const int field_index_r,
                               amrex::MultiFab& mf_t, const int field_index_t);

        void FABZForwardTransform(amrex::MFIter const & mfi,
                                  amrex::MultiFab const & tempHTransformedSplit,
                                  int field_index, const bool is_nodal_z);
        void FABZBackwardTransform (amrex::MFIter const & mfi, const int field_index,
                                    amrex::MultiFab & tempHTransformedSplit,
                                    const bool is_nodal_z);

        // Returns an array that holds the kr for all of the modes
        HankelTransform::RealVector const & getKrArray(amrex::MFIter const & mfi) const {return hankeltransformer[mfi].getKrArray();}

        // `fields` stores fields in spectral space, as multicomponent FabArray
        SpectralField fields;

        int n_rz_azimuthal_modes;

    private:

        // tmpRealField and tmpSpectralField store fields
        // right before/after the Fourier transform
        SpectralField tempHTransformed; // contains Complexes
        SpectralField tmpSpectralField; // contains Complexes
        FFTplans forward_plan, backward_plan;
        // Correcting "shift" factors when performing FFT from/to
        // a cell-centered grid in real space, instead of a nodal grid
        SpectralShiftFactor zshift_FFTfromCell, zshift_FFTtoCell;
        HankelTransformer hankeltransformer;

};

#endif // WARPX_SPECTRAL_FIELD_DATA_HANKEL_H_
