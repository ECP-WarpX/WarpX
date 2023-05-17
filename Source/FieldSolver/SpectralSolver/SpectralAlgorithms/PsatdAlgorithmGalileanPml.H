/* Copyright 2019-2023 Axel Huebl, Remi Lehe, Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_GALILEAN_PML_H_
#define WARPX_PSATD_ALGORITHM_GALILEAN_PML_H_

#include "SpectralBaseAlgorithm.H"

#if WARPX_USE_PSATD

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PsatdAlgorithmGalileanPml : public SpectralBaseAlgorithm
{
    public:

        /**
         * \brief Constructor of the class PsatdAlgorithmGalileanPml
         *
         * \param[in] spectral_kspace spectral space
         * \param[in] dm distribution mapping
         * \param[in] spectral_index object containing indices to access data in spectral space
         * \param[in] norder_x order of the spectral solver along x
         * \param[in] norder_y order of the spectral solver along y
         * \param[in] norder_z order of the spectral solver along z
         * \param[in] grid_type type of grid (collocated or not)
         * \param[in] v_galilean Galilean velocity (three-dimensional array)
         * \param[in] dt time step of the simulation
         */
        PsatdAlgorithmGalileanPml (
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const SpectralFieldIndex& spectral_index,
            const int norder_x,
            const int norder_y,
            const int norder_z,
            const short grid_type,
            const amrex::Vector<amrex::Real>& v_galilean,
            const amrex::Real dt);

        /**
         * \brief Initializes the coefficients used in \c pushSpectralFields to update the E and B fields
         *
         * \param[in] spectral_kspace spectral space
         * \param[in] dm distribution mapping
         */
        void InitializeSpectralCoefficients (
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm);

        /**
         * \brief Updates the E and B fields in spectral space, according to the relevant PSATD equations
         *
         * \param[in,out] f all the fields in spectral space
         */
        virtual void pushSpectralFields (SpectralFieldData& f) const override final;

        /**
         * \brief Virtual function for current correction in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c CurrentCorrection in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         */
        virtual void CurrentCorrection (SpectralFieldData& field_data) override final;

        /**
         * \brief Virtual function for Vay current deposition in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c VayDeposition in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         */
        virtual void VayDeposition (SpectralFieldData& field_data) override final;

    private:

        SpectralRealCoefficients C_coef, S_ck_coef;
        SpectralComplexCoefficients T2_coef;

        // Centered modified finite-order k vectors
        KVectorComponent modified_kx_vec_centered;
#if defined(WARPX_DIM_3D)
        KVectorComponent modified_ky_vec_centered;
#endif
        KVectorComponent modified_kz_vec_centered;

        // Other private member variables
        amrex::Vector<amrex::Real> m_v_galilean;
        amrex::Real m_dt;
};

#endif // WARPX_USE_PSATD
#endif // WARPX_PSATD_ALGORITHM_GALILEAN_PML_H_
