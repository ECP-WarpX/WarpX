/* Copyright 2019 Remi Lehe, Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_BASE_ALGORITHM_H_
#define WARPX_SPECTRAL_BASE_ALGORITHM_H_

#include "FieldSolver/SpectralSolver/SpectralKSpace.H"
#include "FieldSolver/SpectralSolver/SpectralFieldData.H"

#include <AMReX.H>

#if WARPX_USE_PSATD

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 *
 * `SpectralBaseAlgorithm` is only a base class and cannot be used directly.
 * Instead use its subclasses, which implement the specific field update
 * equations for a given spectral algorithm.
 */
class SpectralBaseAlgorithm
{
    public:
        // Virtual member function ; meant to be overridden in subclasses
        virtual void pushSpectralFields(SpectralFieldData& f) const = 0;
        virtual int getRequiredNumberOfFields() const = 0;
        // The destructor should also be a virtual function, so that
        // a pointer to subclass of `SpectraBaseAlgorithm` actually
        // calls the subclass's destructor.
        virtual ~SpectralBaseAlgorithm() {}

        /**
         * \brief Virtual function for current correction in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This virtual function is pure and must be defined in derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         * \param[in]     rho        Unique pointer to \c MultiFab storing the charge density
         */
        virtual void CurrentCorrection (const int lev,
                                        SpectralFieldData& field_data,
                                        std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                        const std::unique_ptr<amrex::MultiFab>& rho) = 0;

        /**
         * \brief Virtual function for Vay current deposition in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This virtual function is pure and must be defined in derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         */
        virtual void VayDeposition (const int lev,
                                    SpectralFieldData& field_data,
                                    std::array<std::unique_ptr<amrex::MultiFab>,3>& current) = 0;

        /**
         * \brief Compute spectral divergence of E
         */
        void ComputeSpectralDivE ( const int lev,
                                   SpectralFieldData& field_data,
                                   const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                                   amrex::MultiFab& divE );

    protected: // Meant to be used in the subclasses

        using SpectralRealCoefficients = \
            amrex::FabArray< amrex::BaseFab <amrex::Real> >;
        using SpectralComplexCoefficients = \
            amrex::FabArray< amrex::BaseFab <Complex> >;

        // Constructor
        SpectralBaseAlgorithm(const SpectralKSpace& spectral_kspace,
                              const amrex::DistributionMapping& dm,
                              const int norder_x, const int norder_y,
                              const int norder_z, const bool nodal)
          // Compute and assign the modified k vectors
          : modified_kx_vec(spectral_kspace.getModifiedKComponent(dm,0,norder_x,nodal)),
#if (AMREX_SPACEDIM==3)
            modified_ky_vec(spectral_kspace.getModifiedKComponent(dm,1,norder_y,nodal)),
            modified_kz_vec(spectral_kspace.getModifiedKComponent(dm,2,norder_z,nodal))
#else
            modified_kz_vec(spectral_kspace.getModifiedKComponent(dm,1,norder_z,nodal))
#endif
          {
#if (AMREX_SPACEDIM!=3)
            amrex::ignore_unused(norder_y);
#endif
          }

        // Modified finite-order vectors
        KVectorComponent modified_kx_vec;
#if (AMREX_SPACEDIM==3)
        KVectorComponent modified_ky_vec;
#endif
        KVectorComponent modified_kz_vec;
};

#endif // WARPX_USE_PSATD
#endif // WARPX_SPECTRAL_BASE_ALGORITHM_H_
