/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_FIELD_DATA_RZ_H_
#define WARPX_SPECTRAL_FIELD_DATA_RZ_H_

#include "SpectralKSpaceRZ.H"
#include "SpectralFieldData.H"
#include "SpectralHankelTransform/SpectralHankelTransformer.H"
#include "SpectralBinomialFilter.H"
#include <AMReX_MultiFab.H>

/* \brief Class that stores the fields in spectral space, and performs the
 *  Fourier transforms between real space and spectral space
 */
class SpectralFieldDataRZ
{

    public:

        // Define the FFTplans type, which holds one fft plan per box
        // (plans are only initialized for the boxes that are owned by
        // the local MPI rank)
#if defined(AMREX_USE_CUDA)
        using FFTplans = amrex::LayoutData<cufftHandle>;
#elif defined(AMREX_USE_HIP)
        using FFTplans = amrex::LayoutData<rocfft_plan>;
#else
        using FFTplans = amrex::LayoutData<fftw_plan>;
#endif
        // Similarly, define the Hankel transformers and filter for each box.
        using MultiSpectralHankelTransformer = amrex::LayoutData<SpectralHankelTransformer>;

        using BinomialFilter = amrex::LayoutData<SpectralBinomialFilter>;

        SpectralFieldDataRZ (const amrex::BoxArray& realspace_ba,
                             const SpectralKSpaceRZ& k_space,
                             const amrex::DistributionMapping& dm,
                             const int n_field_required,
                             const int n_modes,
                             const int lev);
        SpectralFieldDataRZ () = default; // Default constructor
        SpectralFieldDataRZ& operator=(SpectralFieldDataRZ&& field_data) = default;
        ~SpectralFieldDataRZ ();

        void ForwardTransform (const amrex::MultiFab& mf, const int field_index,
                               const int i_comp=0);
        void ForwardTransform (const amrex::MultiFab& mf_r, const int field_index_r,
                               const amrex::MultiFab& mf_t, const int field_index_t);
        void BackwardTransform (amrex::MultiFab& mf, const int field_index,
                                const int i_comp=0);
        void BackwardTransform (amrex::MultiFab& mf_r, const int field_index_r,
                                amrex::MultiFab& mf_t, const int field_index_t);

        void FABZForwardTransform (amrex::MFIter const & mfi, amrex::Box const & realspace_bx,
                                   amrex::MultiFab const & tempHTransformedSplit,
                                   int field_index, const bool is_nodal_z);
        void FABZBackwardTransform (amrex::MFIter const & mfi, amrex::Box const & realspace_bx,
                                    const int field_index,
                                    amrex::MultiFab & tempHTransformedSplit,
                                    const bool is_nodal_z);

        void InitFilter (amrex::IntVect const & filter_npass_each_dir, bool const compensation,
                         SpectralKSpaceRZ const & k_space);

        void ApplyFilter (int const field_index);
        void ApplyFilter (int const field_index1, int const field_index2, int const field_index3);

        // Returns an array that holds the kr for all of the modes
        HankelTransform::RealVector const & getKrArray (amrex::MFIter const & mfi) const {
            return multi_spectral_hankel_transformer[mfi].getKrArray();
        }

        // `fields` stores fields in spectral space, as multicomponent FabArray
        SpectralField fields;

        int n_rz_azimuthal_modes;

    private:

        // tempHTransformed and tmpSpectralField store fields
        // right before/after the z Fourier transform
        SpectralField tempHTransformed; // contains Complexes
        SpectralField tmpSpectralField; // contains Complexes
        FFTplans forward_plan, backward_plan;
        // Correcting "shift" factors when performing FFT from/to
        // a cell-centered grid in real space, instead of a nodal grid
        SpectralShiftFactor zshift_FFTfromCell, zshift_FFTtoCell;
        MultiSpectralHankelTransformer multi_spectral_hankel_transformer;
        BinomialFilter binomialfilter;

};

#endif // WARPX_SPECTRAL_FIELD_DATA_RZ_H_
