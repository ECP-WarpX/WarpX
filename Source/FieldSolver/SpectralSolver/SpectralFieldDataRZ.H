/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_FIELD_DATA_RZ_H_
#define WARPX_SPECTRAL_FIELD_DATA_RZ_H_

#include "SpectralBinomialFilter.H"
#include "SpectralFieldData.H"
#include "SpectralHankelTransform/SpectralHankelTransformer.H"
#include "SpectralKSpaceRZ.H"

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

        SpectralFieldDataRZ (const int lev,
                             const amrex::BoxArray& realspace_ba,
                             const SpectralKSpaceRZ& k_space,
                             const amrex::DistributionMapping& dm,
                             const int n_field_required,
                             const int n_modes);
        SpectralFieldDataRZ () = default; // Default constructor
        SpectralFieldDataRZ& operator=(SpectralFieldDataRZ&& field_data) = default;
        ~SpectralFieldDataRZ ();

        void ForwardTransform (const int lev, const amrex::MultiFab& mf, const int field_index,
                               const int i_comp=0);
        void ForwardTransform (const int lev, const amrex::MultiFab& mf_r, const int field_index_r,
                               const amrex::MultiFab& mf_t, const int field_index_t);
        void BackwardTransform (const int lev, amrex::MultiFab& mf, const int field_index,
                                const int i_comp=0);
        void BackwardTransform (const int lev, amrex::MultiFab& mf_r, const int field_index_r,
                                amrex::MultiFab& mf_t, const int field_index_t);

        void FABZForwardTransform (amrex::MFIter const & mfi, amrex::Box const & realspace_bx,
                                   amrex::MultiFab const & tempHTransformedSplit,
                                   int field_index, const bool is_nodal_z);
        void FABZBackwardTransform (amrex::MFIter const & mfi, amrex::Box const & realspace_bx,
                                    const int field_index,
                                    amrex::MultiFab & tempHTransformedSplit,
                                    const bool is_nodal_z);

        /**
         * \brief Copy spectral data from component \c src_comp to component \c dest_comp
         *        of \c fields
         *
         * \param[in] src_comp  component of the source FabArray from which the data are copied
         * \param[in] dest_comp component of the destination FabArray where the data are copied
         */
        void CopySpectralDataComp (const int src_comp, const int dest_comp)
        {
            // The last two arguments represent the number of components and
            // the number of ghost cells to perform this operation
            Copy(fields, fields, src_comp, dest_comp, n_rz_azimuthal_modes, 0);
        }

        /**
         * \brief Set to zero the data on component \c icomp of \c fields
         *
         * \param[in] icomp component of the FabArray where the data are set to zero
         */
        void ZeroOutDataComp(const int icomp)
        {
            // The last argument represents the number of components to perform this operation
            fields.setVal(0., icomp, n_rz_azimuthal_modes);
        }

        /**
         * \brief Scale the data on component \c icomp of \c fields by a given scale factor
         *
         * \param[in] icomp component of the FabArray where the data are scaled
         * \param[in] scale_factor scale factor to use for scaling
         */
        void ScaleDataComp(const int icomp, const amrex::Real scale_factor)
        {
            // The last argument represents the number of components to perform this operation
            fields.mult(scale_factor, icomp, n_rz_azimuthal_modes);
        }

        void InitFilter (amrex::IntVect const & filter_npass_each_dir, bool const compensation,
                         SpectralKSpaceRZ const & k_space);

        void ApplyFilter (const int lev, int const field_index);
        void ApplyFilter (const int lev, int const field_index1,
                          int const field_index2, int const field_index3);

        // Returns an array that holds the kr for all of the modes
        HankelTransform::RealVector const & getKrArray (amrex::MFIter const & mfi) const {
            return multi_spectral_hankel_transformer[mfi].getKrArray();
        }

        //! `fields` stores fields in spectral space, as multicomponent FabArray
        SpectralField fields;

        //! Number of modes for the RZ multi-mode version, see WarpX::n_rz_azimuthal_modes
        int n_rz_azimuthal_modes;
        //! Number of MultiFab components, see WarpX::ncomps
        int m_ncomps;

    private:

        SpectralFieldIndex m_spectral_index;
        int m_n_fields;

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
