/* Copyright 2019
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_H_
#define WARPX_PSATD_ALGORITHM_H_

#include "FieldSolver/SpectralSolver/SpectralFieldData.H"
#include "FieldSolver/SpectralSolver/SpectralKSpace.H"
#include "SpectralBaseAlgorithm.H"

#include <AMReX_Array.H>
#include <AMReX_Config.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

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
         * \param[in] spectral_index object containing indices to access data in spectral space
         * \param[in] norder_x order of the spectral solver along x
         * \param[in] norder_y order of the spectral solver along y
         * \param[in] norder_z order of the spectral solver along z
         * \param[in] nodal whether the E and B fields are defined on a fully nodal grid or a Yee grid
         * \param[in] v_galilean Galilean velocity (three-dimensional array)
         * \param[in] dt time step of the simulation
         * \param[in] update_with_rho whether the update equation for E uses rho or not
         * \param[in] time_averaging whether to use time averaging for large time steps
         * \param[in] J_linear_in_time whether to use two currents computed at the beginning and the end
         *            of the time interval (instead of using one current computed at half time)
         */
        PsatdAlgorithm (
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const SpectralFieldIndex& spectral_index,
            const int norder_x,
            const int norder_y,
            const int norder_z,
            const bool nodal,
            const amrex::IntVect& fill_guards,
            const amrex::Vector<amrex::Real>& v_galilean,
            const amrex::Real dt,
            const bool update_with_rho,
            const bool time_averaging,
            const bool J_linear_in_time,
            const bool dive_cleaning,
            const bool divb_cleaning);

        /**
         * \brief Updates the E and B fields in spectral space, according to the relevant PSATD equations
         *
         * \param[in,out] f all the fields in spectral space
         */
        virtual void pushSpectralFields (SpectralFieldData& f) const override final;

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
         * \brief Initialize additional coefficients used in \c pushSpectralFields to update E,B,
         *        required only when using time averaging with the assumption that J is linear in time
         *
         * \param[in] spectral_kspace spectral space
         * \param[in] dm distribution mapping
         * \param[in] dt time step of the simulation
         */
        void InitializeSpectralCoefficientsAvgLin (
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

        SpectralComplexCoefficients X5_coef, X6_coef;

        // These real and complex coefficients are allocated only with averaged Galilean PSATD
        SpectralComplexCoefficients Psi1_coef, Psi2_coef, Y1_coef, Y2_coef, Y3_coef, Y4_coef;

        SpectralFieldIndex m_spectral_index;

        // Centered modified finite-order k vectors
        KVectorComponent modified_kx_vec_centered;
#if (AMREX_SPACEDIM==3)
        KVectorComponent modified_ky_vec_centered;
#endif
        KVectorComponent modified_kz_vec_centered;

        // Other member variables
        amrex::Vector<amrex::Real> m_v_galilean;
        amrex::Real m_dt;
        bool m_update_with_rho;
        bool m_time_averaging;
        bool m_J_linear_in_time;
        bool m_dive_cleaning;
        bool m_divb_cleaning;
        bool m_is_galilean;
};
#endif // WARPX_USE_PSATD
#endif // WARPX_PSATD_ALGORITHM_H_
