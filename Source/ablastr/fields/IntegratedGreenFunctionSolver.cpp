/* Copyright 2019-2024 Arianna Formenti, Remi Lehe
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "IntegratedGreenFunctionSolver.H"

#include <ablastr/constant.H>
#include <ablastr/warn_manager/WarnManager.H>
#include <ablastr/math/fft/AnyFFT.H>

#include <AMReX_Array4.H>
#include <AMReX_BaseFab.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <array>


namespace ablastr::fields {

void
computePhiIGF ( amrex::MultiFab const & rho,
                amrex::MultiFab & phi,
                std::array<amrex::Real, 3> const & cell_size,
                amrex::BoxArray const & ba )
{
    using namespace amrex::literals;

    // Define box that encompasses the full domain
    amrex::Box domain = ba.minimalBox();
    domain.surroundingNodes(); // get nodal points, since `phi` and `rho` are nodal
    domain.grow( phi.nGrowVect() ); // include guard cells

    int const nx = domain.length(0);
    int const ny = domain.length(1);
    int const nz = domain.length(2);

    // Allocate 2x wider arrays for the convolution of rho with the Green function
    // This also defines the box arrays for the global FFT: contains only one box;
    amrex::Box const realspace_box = amrex::Box(
        {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)},
        {2*nx-1+domain.smallEnd(0), 2*ny-1+domain.smallEnd(1), 2*nz-1+domain.smallEnd(2)},
        amrex::IntVect::TheNodeVector() );
    amrex::BoxArray const realspace_ba = amrex::BoxArray( realspace_box );
    amrex::Box const spectralspace_box = amrex::Box(
        {0,0,0},
        {nx, 2*ny-1, 2*nz-1},
        amrex::IntVect::TheNodeVector() );
    amrex::BoxArray const spectralspace_ba = amrex::BoxArray( spectralspace_box );
    // Define a distribution mapping for the global FFT, with only one box
    amrex::DistributionMapping dm_global_fft;
    dm_global_fft.define( realspace_ba );
    // Allocate required arrays
    amrex::MultiFab tmp_rho = amrex::MultiFab(realspace_ba, dm_global_fft, 1, 0);
    tmp_rho.setVal(0);
    amrex::MultiFab tmp_G = amrex::MultiFab(realspace_ba, dm_global_fft, 1, 0);
    tmp_G.setVal(0);
    // Allocate corresponding arrays in Fourier space
    using SpectralField = amrex::FabArray< amrex::BaseFab< amrex::GpuComplex< amrex::Real > > >;
    SpectralField tmp_rho_fft = SpectralField( spectralspace_ba, dm_global_fft, 1, 0 );
    SpectralField tmp_G_fft = SpectralField( spectralspace_ba, dm_global_fft, 1, 0 );

    // Copy from rho to tmp_rho
    tmp_rho.ParallelCopy( rho, 0, 0, 1, amrex::IntVect::TheZeroVector(), amrex::IntVect::TheZeroVector() );

    // Compute the integrated Green function
    {
    BL_PROFILE("Initialize Green function");
    amrex::BoxArray const domain_ba = amrex::BoxArray( domain );
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(domain_ba, dm_global_fft,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::Box const bx = mfi.tilebox();

        amrex::IntVect const lo = realspace_box.smallEnd();
        amrex::IntVect const hi = realspace_box.bigEnd();

        // Fill values of the Green function
        amrex::Real const dx = cell_size[0];
        amrex::Real const dy = cell_size[1];
        amrex::Real const dz = cell_size[2];
        amrex::Array4<amrex::Real> const tmp_G_arr = tmp_G.array(mfi);
        amrex::ParallelFor( bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                int const i0 = i - lo[0];
                int const j0 = j - lo[1];
                int const k0 = k - lo[2];
                amrex::Real const x = i0*dx;
                amrex::Real const y = j0*dy;
                amrex::Real const z = k0*dz;

                amrex::Real const G_value = 1._rt/(4._rt*ablastr::constant::math::pi*ablastr::constant::SI::ep0) * (
                    IntegratedPotential( x+0.5_rt*dx, y+0.5_rt*dy, z+0.5_rt*dz )
                  - IntegratedPotential( x-0.5_rt*dx, y+0.5_rt*dy, z+0.5_rt*dz )
                  - IntegratedPotential( x+0.5_rt*dx, y-0.5_rt*dy, z+0.5_rt*dz )
                  - IntegratedPotential( x+0.5_rt*dx, y+0.5_rt*dy, z-0.5_rt*dz )
                  + IntegratedPotential( x+0.5_rt*dx, y-0.5_rt*dy, z-0.5_rt*dz )
                  + IntegratedPotential( x-0.5_rt*dx, y+0.5_rt*dy, z-0.5_rt*dz )
                  + IntegratedPotential( x-0.5_rt*dx, y-0.5_rt*dy, z+0.5_rt*dz )
                  - IntegratedPotential( x-0.5_rt*dx, y-0.5_rt*dy, z-0.5_rt*dz )
                );

                tmp_G_arr(i,j,k) = G_value;
                // Fill the rest of the array by periodicity
                if (i0>0) {tmp_G_arr(hi[0]+1-i0, j         , k         ) = G_value;}
                if (j0>0) {tmp_G_arr(i         , hi[1]+1-j0, k         ) = G_value;}
                if (k0>0) {tmp_G_arr(i         , j         , hi[2]+1-k0) = G_value;}
                if ((i0>0)&&(j0>0)) {tmp_G_arr(hi[0]+1-i0, hi[1]+1-j0, k         ) = G_value;}
                if ((j0>0)&&(k0>0)) {tmp_G_arr(i         , hi[1]+1-j0, hi[2]+1-k0) = G_value;}
                if ((i0>0)&&(k0>0)) {tmp_G_arr(hi[0]+1-i0, j         , hi[2]+1-k0) = G_value;}
                if ((i0>0)&&(j0>0)&&(k0>0)) {tmp_G_arr(hi[0]+1-i0, hi[1]+1-j0, hi[2]+1-k0) = G_value;}
            }
        );
    }
    }

    // Perform forward FFTs
    auto forward_plan_rho = ablastr::math::anyfft::FFTplans(spectralspace_ba, dm_global_fft);
    auto forward_plan_G = ablastr::math::anyfft::FFTplans(spectralspace_ba, dm_global_fft);
    // Loop over boxes perform FFTs
    for ( amrex::MFIter mfi(realspace_ba, dm_global_fft); mfi.isValid(); ++mfi ){

        // Note: the size of the real-space box and spectral-space box
        // differ when using real-to-complex FFT. When initializing
        // the FFT plan, the valid dimensions are those of the real-space box.
        const amrex::IntVect fft_size = realspace_ba[mfi].length();

        // FFT of rho
        forward_plan_rho[mfi] = ablastr::math::anyfft::CreatePlan(
            fft_size, tmp_rho[mfi].dataPtr(),
            reinterpret_cast<ablastr::math::anyfft::Complex*>(tmp_rho_fft[mfi].dataPtr()),
            ablastr::math::anyfft::direction::R2C, AMREX_SPACEDIM);
        ablastr::math::anyfft::Execute(forward_plan_rho[mfi]);

        // FFT of G
        forward_plan_G[mfi] = ablastr::math::anyfft::CreatePlan(
            fft_size, tmp_G[mfi].dataPtr(),
            reinterpret_cast<ablastr::math::anyfft::Complex*>(tmp_G_fft[mfi].dataPtr()),
            ablastr::math::anyfft::direction::R2C, AMREX_SPACEDIM);
        ablastr::math::anyfft::Execute(forward_plan_G[mfi]);

    }

    // Multiply tmp_G_fft and tmp_rho_fft in spectral space
    // Store the result in-place in Gtmp_G_fft, to save memory
    amrex::Multiply( tmp_G_fft, tmp_rho_fft, 0, 0, 1, 0);

    // Perform inverse FFT
    auto backward_plan = ablastr::math::anyfft::FFTplans(spectralspace_ba, dm_global_fft);
    // Loop over boxes perform FFTs
    for ( amrex::MFIter mfi(spectralspace_ba, dm_global_fft); mfi.isValid(); ++mfi ){

        // Note: the size of the real-space box and spectral-space box
        // differ when using real-to-complex FFT. When initializing
        // the FFT plan, the valid dimensions are those of the real-space box.
        const amrex::IntVect fft_size = realspace_ba[mfi].length();

        // Inverse FFT: is done in-place, in the array of G
        backward_plan[mfi] = ablastr::math::anyfft::CreatePlan(
            fft_size, tmp_G[mfi].dataPtr(),
            reinterpret_cast<ablastr::math::anyfft::Complex*>( tmp_G_fft[mfi].dataPtr()),
            ablastr::math::anyfft::direction::C2R, AMREX_SPACEDIM);
        ablastr::math::anyfft::Execute(backward_plan[mfi]);
    }
    // Normalize, since (FFT + inverse FFT) results in a factor N
    const amrex::Real normalization = 1._rt / realspace_box.numPts();
    tmp_G.mult( normalization );

    // Copy from tmp_G to phi
    phi.ParallelCopy( tmp_G, 0, 0, 1, amrex::IntVect::TheZeroVector(), phi.nGrowVect() );

    // Loop to destroy FFT plans
    for ( amrex::MFIter mfi(spectralspace_ba, dm_global_fft); mfi.isValid(); ++mfi ){
        ablastr::math::anyfft::DestroyPlan(forward_plan_G[mfi]);
        ablastr::math::anyfft::DestroyPlan(forward_plan_rho[mfi]);
        ablastr::math::anyfft::DestroyPlan(backward_plan[mfi]);
    }
}
} // namespace ablastr::fields
