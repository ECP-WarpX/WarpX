/* Copyright 2019 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_RZ_H_
#define WARPX_PSATD_ALGORITHM_RZ_H_

#include "SpectralBaseAlgorithmRZ.H"

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PsatdAlgorithmRZ : public SpectralBaseAlgorithmRZ
{

    public:
        PsatdAlgorithmRZ(SpectralKSpaceRZ const & spectral_kspace,
                         amrex::DistributionMapping const & dm,
                         int const n_rz_azimuthal_modes, int const norder_z,
                         bool const nodal, amrex::Real const dt_step);
        // Redefine functions from base class
        virtual void pushSpectralFields(SpectralFieldDataRZ & f) override final;
        virtual int getRequiredNumberOfFields() const override final {
            return SpectralFieldIndex::n_fields;
        }

        void InitializeSpectralCoefficients(SpectralFieldDataRZ const & f);

        /**
         * \brief Virtual function for current correction in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c CurrentCorrection in the
         * base class \c SpectralBaseAlgorithmRZ and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         * \param[in]     rho        Unique pointer to \c MultiFab storing the charge density
         */
        virtual void CurrentCorrection (SpectralFieldDataRZ& field_data,
                                        std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                        const std::unique_ptr<amrex::MultiFab>& rho) override final;

        /**
         * \brief Virtual function for Vay current deposition in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c VayDeposition in the
         * base class \c SpectralBaseAlgorithmRZ and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         */
        virtual void VayDeposition (SpectralFieldDataRZ& field_data,
                                    std::array<std::unique_ptr<amrex::MultiFab>,3>& current) override final;

    private:
        bool coefficients_initialized;
        // Note that dt is saved to use in InitializeSpectralCoefficients
        amrex::Real m_dt;
        SpectralRealCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_ALGORITHM_RZ_H_
