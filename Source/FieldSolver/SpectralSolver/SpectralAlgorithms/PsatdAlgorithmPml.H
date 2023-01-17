/* Copyright 2019 Axel Huebl, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_PML_H_
#define WARPX_PSATD_ALGORITHM_PML_H_

#include "SpectralBaseAlgorithm.H"

#include "FieldSolver/SpectralSolver/SpectralFieldData_fwd.H"
#include "FieldSolver/SpectralSolver/SpectralKSpace_fwd.H"

#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

#if WARPX_USE_PSATD

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PsatdAlgorithmPml : public SpectralBaseAlgorithm
{
    public:
        PsatdAlgorithmPml(const SpectralKSpace& spectral_kspace,
                          const amrex::DistributionMapping& dm,
                          const SpectralFieldIndex& spectral_index,
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
        SpectralFieldIndex m_spectral_index;
        SpectralRealCoefficients C_coef, S_ck_coef, inv_k2_coef;
        amrex::Real m_dt;
        bool m_dive_cleaning;
        bool m_divb_cleaning;
};

#endif // WARPX_USE_PSATD
#endif // WARPX_PSATD_ALGORITHM_PML_H_
