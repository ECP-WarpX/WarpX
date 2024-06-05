/* Copyright 2019-2024 Arianna Formenti, Remi Lehe
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <heffte.h>

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
#include <AMReX_PlotFileUtil.H>

#include <array>
#include <fstream>
#include <iomanip>

namespace ablastr::fields {

void
computePhiIGF ( amrex::MultiFab const & rho,
                amrex::MultiFab & phi,
                std::array<amrex::Real, 3> const & cell_size,
                amrex::BoxArray const & ba)
{
    using namespace amrex::literals;
    int XXX = 1;

    if(XXX){

    // Define box that encompasses the full domain
    amrex::Box domain = ba.minimalBox();
    domain.surroundingNodes(); // get nodal points, since `phi` and `rho` are nodal
    domain.grow( phi.nGrowVect() ); // include guard cells

    int const nx = domain.length(0);
    int const ny = domain.length(1);
    int const nz = domain.length(2);
   
    long npts = domain.numPts();
    amrex::Real sqrtnpts = std::sqrt(npts);

    int const nprocs = amrex::ParallelDescriptor::NProcs();

    // Allocate 2x wider arrays for the convolution of rho with the Green function
    amrex::Box const realspace_box = amrex::Box(
        {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)},
        {2*nx-1+domain.smallEnd(0), 2*ny-1+domain.smallEnd(1), 2*nz-1+domain.smallEnd(2)},
        amrex::IntVect::TheNodeVector());


    amrex::BoxArray realspace_ba;
    amrex::DistributionMapping realspace_dm;
    {
        int realspace_nx = realspace_box.length(0);
        int realspace_ny = realspace_box.length(1);
        int realspace_nz = realspace_box.length(2);
        int minsize_z = realspace_nz / nprocs;
        int nleft_z = realspace_nz - minsize_z*nprocs;


        AMREX_ALWAYS_ASSERT(realspace_nz >= nprocs);
        // We are going to split realspace_box in such a way that the first
        // nleft boxes has minsize_z+1 nodes and the others minsize
        // nodes. We do it this way instead of BoxArray::maxSize to make
        // sure there are exactly nprocs boxes and there are no overlaps.

        amrex::BoxList bl(amrex::IndexType::TheNodeType());
        for (int iproc = 0; iproc < nprocs; ++iproc) {
            int zlo, zhi;
            if (iproc < nleft_z) {
                zlo = iproc*(minsize_z+1);
                zhi = zlo + minsize_z;

            } else {
                zlo = iproc*minsize_z + nleft_z;
                zhi = zlo + minsize_z - 1;

            }
            amrex::Box tbx(amrex::IntVect(0,0,zlo),amrex::IntVect(realspace_nx-1,realspace_ny-1,zhi),amrex::IntVect(1));

            tbx.shift(realspace_box.smallEnd());
            bl.push_back(tbx);


        }

        realspace_ba.define(std::move(bl));
        amrex::Vector<int> pmap(nprocs);
        std::iota(pmap.begin(), pmap.end(), 0);
        realspace_dm.define(std::move(pmap));
}


    // Allocate required arrays
    amrex::MultiFab tmp_rho = amrex::MultiFab(realspace_ba, realspace_dm, 1, 0);
    tmp_rho.setVal(0);
    amrex::MultiFab tmp_G = amrex::MultiFab(realspace_ba, realspace_dm, 1, 0);
    tmp_G.setVal(0);

    // Copy from rho including its ghost cells to tmp_rho
    // w.z.: please check. I think we need to use rho.nGrowVect() otherwise
    // the data outside the ba.minimalBox() will not be copied over.
    tmp_rho.ParallelCopy( rho, 0, 0, 1, rho.nGrowVect(), amrex::IntVect::TheZeroVector() );


    // Compute the integrated Green function
    {
    BL_PROFILE("Initialize Green function");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(tmp_G,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::Box const bx = mfi.tilebox();

        amrex::IntVect const lo = realspace_box.smallEnd();
        amrex::IntVect const hi = realspace_box.bigEnd();

        // Fill values of the Green function
        amrex::Real const dx = cell_size[0];
        amrex::Real const dy = cell_size[1];
        amrex::Real const dz = cell_size[2];

        amrex::Real x_hi = dx*(hi[0]+2);
        amrex::Real y_hi = dy*(hi[1]+2);
        amrex::Real z_hi = dz*(hi[2]+2);

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

                if ((i0< nx)&&(j0< ny)&&(k0< nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x     , y     , z     , dx, dy, dz); }
                if ((i0< nx)&&(j0> ny)&&(k0< nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x     , y_hi-y, z     , dx, dy, dz); }
                if ((i0< nx)&&(j0< ny)&&(k0> nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x     , y     , z_hi-z, dx, dy, dz); }
                if ((i0> nx)&&(j0> ny)&&(k0< nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x_hi-x, y_hi-y, z     , dx, dy, dz); }
                if ((i0< nx)&&(j0> ny)&&(k0> nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x     , y_hi-y, z_hi-z, dx, dy, dz); }
                if ((i0> nx)&&(j0< ny)&&(k0> nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x_hi-x, y     , z_hi-z, dx, dy, dz); }
                if ((i0> nx)&&(j0> ny)&&(k0> nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x_hi-x, y_hi-y, z_hi-z, dx, dy, dz); }
                if ((i0> nx)&&(j0< ny)&&(k0< nz)) { tmp_G_arr(i,j,k) = SumOfIntegratedPotential(x_hi-x, y     , z     , dx, dy, dz); }
         }
      );
    }
<<<<<<< Updated upstream

   /*// Save the MultiFab data to a text file
    std::ofstream outfile("myG.txt");
    outfile << std::setprecision(3); // Set precision for floating-point values

    for (amrex::MFIter mfi(tmp_G, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto& arr = tmp_G.array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    outfile << "i: " << i << " j: " << j << " k: " << k << " G: " << arr(i,j,k) << "\n";
                }
            }
        }
    }
    outfile.close();*/

    amrex::Box const realspace_box1 = amrex::Box(
    {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)},
    {2*nx-1+domain.smallEnd(0), 2*ny-1+domain.smallEnd(1), 2*nz-1+domain.smallEnd(2)},
    amrex::IntVect::TheCellVector());
    amrex::Geometry geom1(realspace_box1);
    amrex::WriteSingleLevelPlotfile("myG", tmp_G, {"G"}, geom1, 0, 0);
=======
>>>>>>> Stashed changes

    }

    // since there is 1 MPI rank per box, here each MPI rank obtains its local box and the associated boxid
    int local_boxid = amrex::ParallelDescriptor::MyProc(); // because of how we made the DistributionMapping
    amrex::Box local_nodal_box = realspace_ba[local_boxid];
    amrex::Box local_box(local_nodal_box.smallEnd(), local_nodal_box.bigEnd());
    local_box.shift(-realspace_box.smallEnd()); // This simplifies the setup because the global lo is zero now

    // Since we the domain decompostion is in the z-direction, setting up
    // c_local_box is simple.
    amrex::Box c_local_box = local_box;
    c_local_box.setBig(0, local_box.length(0)/2+1);

    using SpectralField = amrex::BaseFab< amrex::GpuComplex< amrex::Real > > ;
    SpectralField tmp_rho_fft(c_local_box, 1, amrex::The_Device_Arena());
    SpectralField tmp_G_fft(c_local_box, 1, amrex::The_Device_Arena());
    tmp_rho_fft.shift(realspace_box.smallEnd());
    tmp_G_fft.shift(realspace_box.smallEnd());

#ifdef AMREX_USE_CUDA
    heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
    heffte::fft3d_r2c<heffte::backend::rocfft> fft
#else
    heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1), local_box.smallEnd(2)},
          {local_box.bigEnd(0)  ,local_box.bigEnd(1)  , local_box.bigEnd(2)}},
         {{c_local_box.smallEnd(0),c_local_box.smallEnd(1), c_local_box.smallEnd(2)},
          {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
         0, amrex::ParallelDescriptor::Communicator());

    using heffte_complex = typename heffte::fft_output<amrex::Real>::type;
    heffte_complex* rho_fft_data = (heffte_complex*) tmp_rho_fft.dataPtr();
    heffte_complex* G_fft_data = (heffte_complex*) tmp_G_fft.dataPtr();


    fft.forward(tmp_rho[local_boxid].dataPtr(), rho_fft_data);
    fft.forward(tmp_G[local_boxid].dataPtr(), G_fft_data);

    // PRINT / SAVE THE FFT OF RHO AND/OR G

    // PRINT / SAVE THE FFT OF G
    // **********************************

    amrex::BoxArray fft_ba;
    {
        amrex::BoxList bl(amrex::IndexType::TheNodeType());
        bl.reserve(realspace_ba.size());

        for (int i = 0; i < realspace_ba.size(); ++i) {
            amrex::Box b = realspace_ba[i];

            amrex::Box r_box = b;
            amrex::Box c_box = amrex::coarsen(r_box, amrex::IntVect(AMREX_D_DECL(2,1,1)));

            // this avoids overlap for the cases when one or more r_box's
            // have an even cell index in the hi-x cell
            if (c_box.bigEnd(0) * 2 == r_box.bigEnd(0)) {
                c_box.setBig(0,c_box.bigEnd(0)-1);
            }

            // increase the size of boxes touching the hi-x domain by 1 in x
            // this is an (Nx x Ny x Nz) -> (Nx/2+1 x Ny x Nz) real-to-complex sizing
            if (b.bigEnd(0) == domain.bigEnd(0)) {
                c_box.growHi(0,1);
            }
            bl.push_back(c_box);

        }
        fft_ba.define(std::move(bl));
    }



    // storage for real, imaginary, magnitude, and phase
    amrex::MultiFab fft_data(fft_ba,realspace_dm,4,0);

    // this copies the spectral data into a distributed MultiFab
    for (amrex::MFIter mfi(fft_data); mfi.isValid(); ++mfi) {

        amrex::Array4<amrex::Real> const& data = fft_data.array(mfi);
        amrex::Array4< amrex::GpuComplex<amrex::Real> > spectral = tmp_G_fft.array();

        const amrex::Box& bx = mfi.validbox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            amrex::Real re = spectral(i,j,k).real() / sqrtnpts;
            amrex::Real im = spectral(i,j,k).imag() / sqrtnpts;

            data(i,j,k,0) = re;
            data(i,j,k,1) = im;
            data(i,j,k,2) = std::sqrt(re*re + im*im);

            // Here we want to store the values of the phase angle
            // Avoid division by zero
            if (re == 0.0) {
                if (im == 0.0){
                    data(i,j,k,3) = 0.0;
                } else if (im > 0.0) {
                    data(i,j,k,3) = M_PI/2.0;
                } else {
                    data(i,j,k,3) = -M_PI/2.0;
                }
            } else {
                data(i,j,k,3) = std::atan(im/re);
            }
        });
    }


    // domain for G fft data used to contruct a geometry object
    amrex::Box domain_fft = amrex::coarsen(domain, amrex::IntVect(AMREX_D_DECL(2,1,1)));
    // shrink by 1 in x in case there are an odd number of cells in the x-direction in domain
    if (domain_fft.bigEnd(0) * 2 == domain.bigEnd(0)) {
        domain_fft.setBig(0,domain_fft.bigEnd(0)-1);
    }
    // grow by 1 in the x-direction to match the size of the FFT
    domain_fft.growHi(0,1);

// **********************************

    // Multiply tmp_G_fft and tmp_rho_fft in spectral space
    // Store the result in-place in Gtmp_G_fft, to save memory
    //amrex::Multiply( tmp_G_fft, tmp_rho_fft, 0, 0, 1, 0);
    tmp_G_fft.mult(tmp_rho_fft, 0, 0, 1);

    // PRINT / SAVE G TIMES RHO

    fft.backward(G_fft_data, tmp_G[local_boxid].dataPtr());

     // Normalize, since (FFT + inverse FFT) results in a factor N
    const amrex::Real normalization = 1._rt / realspace_box.numPts();
    tmp_G.mult( normalization );


    // Copy from tmp_G to phi
    phi.ParallelCopy( tmp_G, 0, 0, 1, amrex::IntVect::TheZeroVector(), phi.nGrowVect());


    }
    else{

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
<<<<<<< Updated upstream

    /*// Save the MultiFab data (phi) to a text file
    std::ofstream outfile("remiphi.txt");
    outfile << std::setprecision(3); // Set precision for floating-point values

    for (amrex::MFIter mfi(phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto& arr = phi.array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    outfile << "i: " << i << " j: " << j << " k: " << k << " phi: " << arr(i,j,k) << "\n";
                }
            }
        }
    }
    outfile.close();*/


    amrex::Box const realspace_box1 = amrex::Box(
    {domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)},
    {2*nx-1+domain.smallEnd(0), 2*ny-1+domain.smallEnd(1), 2*nz-1+domain.smallEnd(2)},
    amrex::IntVect::TheCellVector() );
    amrex::Geometry geom1(realspace_box1);
    amrex::WriteSingleLevelPlotfile("remiphi", tmp_G, {"phi"}, geom1, 0, 0);


=======
>>>>>>> Stashed changes

    // Loop to destroy FFT plans
    for ( amrex::MFIter mfi(spectralspace_ba, dm_global_fft); mfi.isValid(); ++mfi ){
        ablastr::math::anyfft::DestroyPlan(forward_plan_G[mfi]);
        ablastr::math::anyfft::DestroyPlan(forward_plan_rho[mfi]);
        ablastr::math::anyfft::DestroyPlan(backward_plan[mfi]);
    }
    }
}
} // namespace ablastr::fields
