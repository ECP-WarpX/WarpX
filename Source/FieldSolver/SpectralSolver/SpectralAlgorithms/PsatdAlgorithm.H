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

        /**
         * \brief Constructor of the class PsatdAlgorithm
         *
         * \param[in] spectral_kspace spectral space
         * \param[in] dm distribution mapping
         * \param[in] norder_x order of the spectral solver along x
         * \param[in] norder_y order of the spectral solver along y
         * \param[in] norder_z order of the spectral solver along z
         * \param[in] nodal whether the E and B fields are defined on a fully nodal grid or a Yee grid
         * \param[in] v_galilean Galilean velocity (three-dimensional array)
         * \param[in] dt time step of the simulation
         * \param[in] update_with_rho whether the update equation for E uses rho or not
         * \param[in] time_averaging whether to use time averaging for large time steps
         */
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

        /**
         * \brief Updates the E and B fields in spectral space, according to the relevant PSATD equations
         *
         * \param[in,out] f all the fields in spectral space
         */
        virtual void pushSpectralFields (SpectralFieldData& f) const override final;

        /**
         * \brief Returns the number of fields stored in spectral space
         */
        virtual int getRequiredNumberOfFields () const override final
        {
            if (m_time_averaging) {
                return SpectralAvgFieldIndex::n_fields;
            } else {
                return SpectralFieldIndex::n_fields;
            }
        }

        /**
         * \brief Initializes the coefficients used in \c pushSpectralFields to update the E and B fields
         *
         * \param[in] spectral_kspace spectral space
         * \param[in] dm distribution mapping
         * \param[in] dt time step of the simulation
         */
        void InitializeSpectralCoefficients (
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const amrex::Real dt);

        /**
         * \brief Initializes additional coefficients used in \c pushSpectralFields to update the E and B fields,
         *        required only when using time averaging with large time steps
         *
         * \param[in] spectral_kspace spectral space
         * \param[in] dm distribution mapping
         * \param[in] dt time step of the simulation
         */
        void InitializeSpectralCoefficientsAveraging (
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
        SpectralComplexCoefficients Psi1_coef, Psi2_coef, Y1_coef, Y2_coef, Y3_coef, Y4_coef;

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
