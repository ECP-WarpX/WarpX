/* Copyright 2019 Maxence Thevenet, Remi Lehe, Revathi Jambunathan, Edoardo Zoni
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_H_
#define WARPX_PSATD_ALGORITHM_H_

#include "SpectralBaseAlgorithm.H"


#if WARPX_USE_PSATD

/**
 * \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PsatdAlgorithm : public SpectralBaseAlgorithm
{

    public:
        PsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                         const amrex::DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal,
                         const amrex::Real dt);
        // Redefine functions from base class
        virtual void pushSpectralFields(SpectralFieldData& f) const override final;
        virtual int getRequiredNumberOfFields() const override final {
            return SpectralFieldIndex::n_fields;
        }

        void InitializeSpectralCoefficients(const SpectralKSpace& spectral_kspace,
                                    const amrex::DistributionMapping& dm,
                                    const amrex::Real dt);

        /**
         * \brief Current correction in Fourier space
         * (equation (19) of https://doi.org/10.1016/j.jcp.2013.03.010).
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Three-dimensional array of unique pointers to
         *                           <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
         *                           (current density)
         * \param[in]     rho        Unique pointer to <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
         *                           (charge density)
         */
        virtual void CurrentCorrection ( SpectralFieldData& field_data,
                                         std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                         std::unique_ptr<amrex::MultiFab>& rho );

    private:
        SpectralRealCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
        amrex::Real m_dt;
};

#endif // WARPX_USE_PSATD
#endif // WARPX_PSATD_ALGORITHM_H_
