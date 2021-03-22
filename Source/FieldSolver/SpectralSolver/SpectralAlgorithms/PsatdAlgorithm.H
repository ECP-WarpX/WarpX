/* Copyright 2019
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_H_
#define WARPX_PSATD_ALGORITHM_H_

#include "SpectralBaseAlgorithm.H"

#if WARPX_USE_PSATD
/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PsatdAlgorithm : public SpectralBaseAlgorithm
{
    public:

        // TODO Add Doxygen docs
        PsatdAlgorithm (
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const int norder_x,
            const int norder_y,
            const int norder_z,
            const bool nodal,
            const amrex::Array<amrex::Real,3>& v_galilean,
            const amrex::Real dt,
            const bool update_with_rho,
            const bool time_averaging);

        // TODO Add Doxygen docs
        virtual void pushSpectralFields (SpectralFieldData& f) const override final;

        // TODO Add Doxygen docs
        virtual int getRequiredNumberOfFields () const override final
        {
            if (m_time_averaging) {
                return SpectralAvgFieldIndex::n_fields;
            } else {
                return SpectralFieldIndex::n_fields;
            }
        }

        // TODO Add Doxygen docs
        void InitializeSpectralCoefficients (
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const amrex::Real dt);

        /**
         * \brief Virtual function for current correction in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c CurrentCorrection in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         * \param[in]     rho        Unique pointer to \c MultiFab storing the charge density
         */
        virtual void CurrentCorrection (
            const int lev,
            SpectralFieldData& field_data,
            std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
            const std::unique_ptr<amrex::MultiFab>& rho) override final;

        /**
         * \brief Virtual function for Vay current deposition in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c VayDeposition in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         */
        virtual void VayDeposition (
            const int lev,
            SpectralFieldData& field_data,
            std::array<std::unique_ptr<amrex::MultiFab>,3>& current) override final;

    private:

        // These real and complex coefficients are always allocated
        SpectralRealCoefficients C_coef, S_ck_coef;
        SpectralComplexCoefficients T2_coef, X1_coef, X2_coef, X3_coef, X4_coef;

        // These real and complex coefficients are allocated only with averaged Galilean PSATD
        SpectralRealCoefficients C1_coef, C3_coef, S1_coef, S3_coef;
        SpectralComplexCoefficients Psi1_coef, Psi2_coef, Psi3_coef, Psi4_coef,
                                    A1_coef, A2_coef, Rhoold_coef, Rhonew_coef, Jcoef_coef;

        // Centered modified finite-order k vectors
        KVectorComponent modified_kx_vec_centered;
#if (AMREX_SPACEDIM==3)
        KVectorComponent modified_ky_vec_centered;
#endif
        KVectorComponent modified_kz_vec_centered;

        // Other member variables
        amrex::Array<amrex::Real,3> m_v_galilean;
        amrex::Real m_dt;
        bool m_update_with_rho;
        bool m_time_averaging;
        bool m_is_galilean;
};
#endif // WARPX_USE_PSATD
#endif // WARPX_PSATD_ALGORITHM_H_
