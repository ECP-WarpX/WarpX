/* Copyright 2019 Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SpectralBaseAlgorithm.H"
#include "FieldSolver/SpectralSolver/SpectralFieldData.H"

#include "Utils/WarpX_Complex.H"

#include <AMReX_BaseFab.H>
#include <AMReX_Config.H>
#include <AMReX_FabArray.H>
#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_GpuLaunchFunctsC.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_PODVector.H>

#include <array>
#include <memory>

using namespace amrex;

/**
 * \brief Constructor
 */
SpectralBaseAlgorithm::SpectralBaseAlgorithm(const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const int norder_x, const int norder_y,
    const int norder_z, const bool nodal):
    // Compute and assign the modified k vectors
        modified_kx_vec(spectral_kspace.getModifiedKComponent(dm,0,norder_x,nodal)),
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

/**
 * \brief Compute spectral divergence of E
 */
void
SpectralBaseAlgorithm::ComputeSpectralDivE (
    const int lev,
    SpectralFieldData& field_data,
    const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
    amrex::MultiFab& divE )
{
    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of E
    field_data.ForwardTransform(lev, *Efield[0], Idx::Ex, 0 );
    field_data.ForwardTransform(lev, *Efield[1], Idx::Ey, 0 );
    field_data.ForwardTransform(lev, *Efield[2], Idx::Ez, 0 );

    // Loop over boxes
    for (MFIter mfi(field_data.fields); mfi.isValid(); ++mfi){

        const Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = field_data.fields[mfi].array();
        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Shortcuts for the components of E
            const Complex Ex = fields(i,j,k,Idx::Ex);
            const Complex Ey = fields(i,j,k,Idx::Ey);
            const Complex Ez = fields(i,j,k,Idx::Ez);
            // k vector values
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            const Complex I = Complex{0,1};

            // div(E) in Fourier space
            fields(i,j,k,Idx::divE) = I*(kx*Ex+ky*Ey+kz*Ez);
        });
    }

    // Backward Fourier transform
    field_data.BackwardTransform(lev, divE, Idx::divE, 0 );
}
