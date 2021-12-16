#ifndef WARPX_COMOVING_PSATD_ALGORITHM_H_
#define WARPX_COMOVING_PSATD_ALGORITHM_H_

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

/* \brief Class that updates the field in spectral space and stores the coefficients
 * of the corresponding update equation, according to the comoving spectral scheme.
 */
class ComovingPsatdAlgorithm : public SpectralBaseAlgorithm
{
    public:

        /**
         * \brief Class constructor
         */
        ComovingPsatdAlgorithm (const SpectralKSpace& spectral_kspace,
                                const amrex::DistributionMapping& dm,
                                const SpectralFieldIndex& spectral_index,
                                const int norder_x,
                                const int norder_y,
                                const int norder_z,
                                const bool nodal,
                                const amrex::IntVect& fill_guards,
                                const amrex::Vector<amrex::Real>& v_comoving,
                                const amrex::Real dt,
                                const bool update_with_rho);

        /**
         * \brief Override the update equations in Fourier space
         */
        virtual void pushSpectralFields (SpectralFieldData& f) const override final;

        /* \brief Initialize the coefficients needed in the update equations
         */
        void InitializeSpectralCoefficients (const SpectralKSpace& spectral_kspace,
                                             const amrex::DistributionMapping& dm,
                                             const amrex::Real dt);

        /**
         * \brief Virtual function for current correction in Fourier space.
         * This function overrides the virtual function \c CurrentCorrection in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in]     lev        The mesh-refinement level
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
         * \brief Virtual function for Vay current deposition in Fourier space.
         * This function overrides the virtual function \c VayDeposition in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in]     lev        The mesh-refinement level
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         */
        virtual void VayDeposition (const int lev,
                                    SpectralFieldData& field_data,
                                    std::array<std::unique_ptr<amrex::MultiFab>,3>& current) override final;

    private:

        // Real and complex spectral coefficients
        SpectralRealCoefficients C_coef, S_ck_coef;
        SpectralComplexCoefficients Theta2_coef, X1_coef, X2_coef, X3_coef, X4_coef;

        SpectralFieldIndex m_spectral_index;

        // k vectors
        KVectorComponent kx_vec;
#if defined(WARPX_DIM_3D)
        KVectorComponent ky_vec;
#endif
        KVectorComponent kz_vec;

        // Additional member variables
        amrex::Vector<amrex::Real> m_v_comoving;
        amrex::Real m_dt;
};

#endif // WARPX_USE_PSATD
#endif // WARPX_COMOVING_PSATD_ALGORITHM_H_
