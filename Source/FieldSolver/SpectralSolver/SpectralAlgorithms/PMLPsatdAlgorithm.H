/* Copyright 2019 Axel Huebl, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PML_PSATD_ALGORITHM_H_
#define WARPX_PML_PSATD_ALGORITHM_H_

#include "SpectralBaseAlgorithm.H"

#if WARPX_USE_PSATD

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PMLPsatdAlgorithm : public SpectralBaseAlgorithm
{
    public:
        PMLPsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                         const amrex::DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal,
                         const amrex::Real dt,
                         const bool dive_cleaning,
                         const bool divb_cleaning);

        void InitializeSpectralCoefficients(
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const amrex::Real dt);

        // Redefine functions from base class
        virtual void pushSpectralFields(SpectralFieldData& f) const override final;
        virtual int getRequiredNumberOfFields() const override final {
            return SpectralPMLIndex::n_fields;
        }

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
        virtual void CurrentCorrection (const int lev,
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
        virtual void VayDeposition (const int lev,
                                    SpectralFieldData& field_data,
                                    std::array<std::unique_ptr<amrex::MultiFab>,3>& current) override final;

    private:
        SpectralRealCoefficients C_coef, S_ck_coef, inv_k2_coef;
        amrex::Real m_dt;
        bool m_dive_cleaning;
        bool m_divb_cleaning;
};

#endif // WARPX_USE_PSATD
#endif // WARPX_PML_PSATD_ALGORITHM_H_
