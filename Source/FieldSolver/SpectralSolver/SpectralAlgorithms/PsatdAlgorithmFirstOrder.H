/* Copyright 2019
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_FIRST_ORDER_H_
#define WARPX_PSATD_ALGORITHM_FIRST_ORDER_H_

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
/*
 * \brief Class that updates the fields in spectral space according to the first-order PSATD equations.
 */
class PsatdAlgorithmFirstOrder : public SpectralBaseAlgorithm
{
    public:

        /**
         * \brief Constructor of the class PsatdAlgorithmFirstOrder
         *
         * \param[in] spectral_kspace spectral space
         * \param[in] dm distribution mapping
         * \param[in] spectral_index object containing indices to access data in spectral space
         * \param[in] norder_x order of the spectral solver along x
         * \param[in] norder_y order of the spectral solver along y
         * \param[in] norder_z order of the spectral solver along z
         * \param[in] nodal whether the E and B fields are defined on a fully nodal grid or a Yee grid
         * \param[in] dt time step of the simulation
         * \param[in] div_cleaning whether to use divergence correction for both E and B (thus, F and G)
         * \param[in] J_in_time time dependency of J (currently supported: constant, linear)
         * \param[in] rho_in_time time dependency of rho (currently supported: constant, linear)
         */
        PsatdAlgorithmFirstOrder (
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const SpectralFieldIndex& spectral_index,
            const int norder_x,
            const int norder_y,
            const int norder_z,
            const bool nodal,
            const amrex::Real dt,
            const bool div_cleaning,
            const int J_in_time,
            const int rho_in_time);

        /**
         * \brief Updates E, B, F, and G fields in spectral space,
         *        according to the first-order PSATD equations
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

        SpectralFieldIndex m_spectral_index;

        // Other member variables
        amrex::Real m_dt;
        bool m_div_cleaning;
        int m_J_in_time;
        int m_rho_in_time;
};
#endif // WARPX_USE_PSATD
#endif // WARPX_PSATD_ALGORITHM_FIRST_ORDER_H_
