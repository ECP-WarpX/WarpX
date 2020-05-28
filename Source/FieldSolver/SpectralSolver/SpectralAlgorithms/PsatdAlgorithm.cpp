/* Copyright 2019 Remi Lehe, Revathi Jambunathan, Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "PsatdAlgorithm.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <cmath>


#if WARPX_USE_PSATD
using namespace amrex;

/**
 * \brief Constructor
 */
PsatdAlgorithm::PsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                         const DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal, const Real dt)
     // Initialize members of base class
     : SpectralBaseAlgorithm( spectral_kspace, dm,
                              norder_x, norder_y, norder_z, nodal )
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X1_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X2_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X3_coef = SpectralRealCoefficients(ba, dm, 1, 0);

    // Initialize coefficients for update equations
    InitializeSpectralCoefficients(spectral_kspace, dm, dt);

    m_dt = dt;
}

/**
 * \brief Advance E and B fields in spectral space (stored in `f`) over one time step
 */
void
PsatdAlgorithm::pushSpectralFields(SpectralFieldData& f) const{

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        const Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();
        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
        Array4<const Real> X1_arr = X1_coef[mfi].array();
        Array4<const Real> X2_arr = X2_coef[mfi].array();
        Array4<const Real> X3_arr = X3_coef[mfi].array();
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
            // Record old values of the fields to be updated
            using Idx = SpectralFieldIndex;
            const Complex Ex_old = fields(i,j,k,Idx::Ex);
            const Complex Ey_old = fields(i,j,k,Idx::Ey);
            const Complex Ez_old = fields(i,j,k,Idx::Ez);
            const Complex Bx_old = fields(i,j,k,Idx::Bx);
            const Complex By_old = fields(i,j,k,Idx::By);
            const Complex Bz_old = fields(i,j,k,Idx::Bz);
            // Shortcut for the values of J and rho
            const Complex Jx = fields(i,j,k,Idx::Jx);
            const Complex Jy = fields(i,j,k,Idx::Jy);
            const Complex Jz = fields(i,j,k,Idx::Jz);
            const Complex rho_old = fields(i,j,k,Idx::rho_old);
            const Complex rho_new = fields(i,j,k,Idx::rho_new);
            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            constexpr Real c2 = PhysConst::c*PhysConst::c;
            constexpr Real inv_ep0 = 1./PhysConst::ep0;
            const Complex I = Complex{0,1};
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Real X1 = X1_arr(i,j,k);
            const Real X2 = X2_arr(i,j,k);
            const Real X3 = X3_arr(i,j,k);


            // Update E (see WarpX online documentation: theory section)
            fields(i,j,k,Idx::Ex) = C*Ex_old
                        + S_ck*(c2*I*(ky*Bz_old - kz*By_old) - inv_ep0*Jx)
                        - I*(X2*rho_new - X3*rho_old)*kx;
            fields(i,j,k,Idx::Ey) = C*Ey_old
                        + S_ck*(c2*I*(kz*Bx_old - kx*Bz_old) - inv_ep0*Jy)
                        - I*(X2*rho_new - X3*rho_old)*ky;
            fields(i,j,k,Idx::Ez) = C*Ez_old
                        + S_ck*(c2*I*(kx*By_old - ky*Bx_old) - inv_ep0*Jz)
                        - I*(X2*rho_new - X3*rho_old)*kz;
            // Update B (see WarpX online documentation: theory section)
            fields(i,j,k,Idx::Bx) = C*Bx_old
                        - S_ck*I*(ky*Ez_old - kz*Ey_old)
                        +   X1*I*(ky*Jz     - kz*Jy);
            fields(i,j,k,Idx::By) = C*By_old
                        - S_ck*I*(kz*Ex_old - kx*Ez_old)
                        +   X1*I*(kz*Jx     - kx*Jz);
            fields(i,j,k,Idx::Bz) = C*Bz_old
                        - S_ck*I*(kx*Ey_old - ky*Ex_old)
                        +   X1*I*(kx*Jy     - ky*Jx);
        });
    }
};

/**
 * \brief Initialize coefficients for update equations
 */
void PsatdAlgorithm::InitializeSpectralCoefficients(const SpectralKSpace& spectral_kspace,
                                    const amrex::DistributionMapping& dm,
                                    const amrex::Real dt)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;
    // Fill them with the right values:
    // Loop over boxes and allocate the corresponding coefficients
    // for each box owned by the local MPI proc
    for (MFIter mfi(ba, dm); mfi.isValid(); ++mfi){

        const Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const Real* modified_kx = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz = modified_kz_vec[mfi].dataPtr();
        // Extract arrays for the coefficients
        Array4<Real> C = C_coef[mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();
        Array4<Real> X1 = X1_coef[mfi].array();
        Array4<Real> X2 = X2_coef[mfi].array();
        Array4<Real> X3 = X3_coef[mfi].array();

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of vector
            const Real k_norm = std::sqrt(
                std::pow(modified_kx[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(modified_ky[j], 2) +
                std::pow(modified_kz[k], 2));
#else
                std::pow(modified_kz[j], 2));
#endif


            // Calculate coefficients
            constexpr Real c = PhysConst::c;
            constexpr Real ep0 = PhysConst::ep0;
            if (k_norm != 0){
                C(i,j,k) = std::cos(c*k_norm*dt);
                S_ck(i,j,k) = std::sin(c*k_norm*dt)/(c*k_norm);
                X1(i,j,k) = (1. - C(i,j,k))/(ep0 * c*c * k_norm*k_norm);
                X2(i,j,k) = (1. - S_ck(i,j,k)/dt)/(ep0 * k_norm*k_norm);
                X3(i,j,k) = (C(i,j,k) - S_ck(i,j,k)/dt)/(ep0 * k_norm*k_norm);
            } else { // Handle k_norm = 0, by using the analytical limit
                C(i,j,k) = 1.;
                S_ck(i,j,k) = dt;
                X1(i,j,k) = 0.5 * dt*dt / ep0;
                X2(i,j,k) = c*c * dt*dt / (6.*ep0);
                X3(i,j,k) = - c*c * dt*dt / (3.*ep0);
            }
        });
     }
}

void
PsatdAlgorithm::CurrentCorrection( SpectralFieldData& field_data,
                                   std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                   const std::unique_ptr<amrex::MultiFab>& rho ) {
    // Profiling
    WARPX_PROFILE( "PsatdAlgorithm::CurrentCorrection" );

    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of J and rho
    field_data.ForwardTransform( *current[0], Idx::Jx, 0 );
    field_data.ForwardTransform( *current[1], Idx::Jy, 0 );
    field_data.ForwardTransform( *current[2], Idx::Jz, 0 );
    field_data.ForwardTransform( *rho, Idx::rho_old, 0 );
    field_data.ForwardTransform( *rho, Idx::rho_new, 1 );

    // Loop over boxes
    for (MFIter mfi(field_data.fields); mfi.isValid(); ++mfi){

        const Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = field_data.fields[mfi].array();
        // Extract pointers for the k vectors
        const Real* const modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* const modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* const modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Local copy of member variables before GPU loop
        const Real dt = m_dt;

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Record old values of the fields to be updated
            using Idx = SpectralFieldIndex;
            // Shortcuts for the values of J and rho
            const Complex Jx = fields(i,j,k,Idx::Jx);
            const Complex Jy = fields(i,j,k,Idx::Jy);
            const Complex Jz = fields(i,j,k,Idx::Jz);
            const Complex rho_old = fields(i,j,k,Idx::rho_old);
            const Complex rho_new = fields(i,j,k,Idx::rho_new);
            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            constexpr Complex I = Complex{0,1};

            const Real k_norm = std::sqrt( kx*kx + ky*ky + kz*kz );

            // div(J) in Fourier space
            const Complex k_dot_J = kx*Jx + ky*Jy + kz*Jz;

            // Correct J
            if ( k_norm != 0 )
            {
                fields(i,j,k,Idx::Jx) = Jx - (k_dot_J-I*(rho_new-rho_old)/dt)*kx/(k_norm*k_norm);
                fields(i,j,k,Idx::Jy) = Jy - (k_dot_J-I*(rho_new-rho_old)/dt)*ky/(k_norm*k_norm);
                fields(i,j,k,Idx::Jz) = Jz - (k_dot_J-I*(rho_new-rho_old)/dt)*kz/(k_norm*k_norm);
            }
        });
    }

    // Backward Fourier transform of J
    field_data.BackwardTransform( *current[0], Idx::Jx, 0 );
    field_data.BackwardTransform( *current[1], Idx::Jy, 0 );
    field_data.BackwardTransform( *current[2], Idx::Jz, 0 );
}

void
PsatdAlgorithm::VayDeposition( SpectralFieldData& field_data,
                               std::array<std::unique_ptr<amrex::MultiFab>,3>& current ) {
    // Profiling
    WARPX_PROFILE( "PsatdAlgorithm::VayDeposition" );

    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of D (temporarily stored in current)
    // (units multiplied by [L])
    field_data.ForwardTransform( *current[0], Idx::Jx, 0 );
    field_data.ForwardTransform( *current[1], Idx::Jy, 0 );
    field_data.ForwardTransform( *current[2], Idx::Jz, 0 );

    // Index type of D and J (either nodal or cell-centered)
    const amrex::IntVect stag_jx = current[0]->ixType().toIntVect();
    const amrex::IntVect stag_jy = current[1]->ixType().toIntVect();
    const amrex::IntVect stag_jz = current[2]->ixType().toIntVect();

    // Index of z direction in 2D and 3D
    constexpr int zdir = AMREX_SPACEDIM-1;

    const amrex::Real dx = m_dx[0];
#if (AMREX_SPACEDIM==3)
    const amrex::Real dy = m_dx[1];
#endif
    const amrex::Real dz = m_dx[zdir];

    // Loop over boxes
    for (amrex::MFIter mfi(field_data.fields); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> fields = field_data.fields[mfi].array();

        // Extract pointers for the k vectors
        const amrex::Real* const kx_arr = kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const amrex::Real* const ky_arr = ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* const kz_arr = kz_vec[mfi].dataPtr();

        // Extract pointers for the modified k vectors
        const amrex::Real* const modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const amrex::Real* const modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* const modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            // Record old values of the fields to be updated
            using Idx = SpectralFieldIndex;

            // Shortcuts for the values of J and rho
            const Complex Dx = fields(i,j,k,Idx::Jx);
            const Complex Dy = fields(i,j,k,Idx::Jy);
            const Complex Dz = fields(i,j,k,Idx::Jz);

            // k vector values
            const amrex::Real kx = kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const amrex::Real ky = ky_arr[j];
            const amrex::Real kz = kz_arr[k];
#else
            constexpr amrex::Real ky = 0;
            const     amrex::Real kz = kz_arr[j];
#endif
            // Imaginary unit
            constexpr Complex I = Complex{0,1};

            // Modified k vector values
            const amrex::Real kx_mod = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const amrex::Real ky_mod = modified_ky_arr[j];
            const amrex::Real kz_mod = modified_kz_arr[k];
#else
            constexpr amrex::Real ky_mod = 0;
            const     amrex::Real kz_mod = modified_kz_arr[j];
#endif

            // Compute Jx cell-centered or nodal
            // (units multiplied by [L] again)
            if ( kx_mod != 0.0_rt ) {
                if      ( stag_jx[0] == 0 ) fields(i,j,k,Idx::Jx) = I*Dx/kx_mod*exp(I*kx*dx*0.5_rt);
                else if ( stag_jx[0] == 1 ) fields(i,j,k,Idx::Jx) = I*Dx/kx_mod;
            }

#if (AMREX_SPACEDIM==3)
            // Compute Jy cell-centered or nodal
            // (units multiplied by [L] again)
            if ( ky_mod != 0.0_rt ) {
                if      ( stag_jy[1] == 0 ) fields(i,j,k,Idx::Jy) = I*Dy/ky_mod*exp(I*ky*dy*0.5_rt);
                else if ( stag_jy[1] == 1 ) fields(i,j,k,Idx::Jy) = I*Dy/ky_mod;
            }
#endif

            // Compute Jz cell-centered or nodal
            // (units multiplied by [L] again)
            if ( kz_mod != 0.0_rt ) {
                if       ( stag_jz[zdir] == 0 ) fields(i,j,k,Idx::Jz) = I*Dz/kz_mod*exp(I*kz*dz*0.5_rt);
                else  if ( stag_jz[zdir] == 1 ) fields(i,j,k,Idx::Jz) = I*Dz/kz_mod;
            }
        } );
    }

    // Alias MultiFab for D before overwriting
    amrex::MultiFab Dx_mf( *current[0], amrex::make_alias, 0, 1 );
    amrex::MultiFab Dy_mf( *current[1], amrex::make_alias, 0, 1 );
    amrex::MultiFab Dz_mf( *current[2], amrex::make_alias, 0, 1 );

    // Temporary MultiFabs for cumulative sums
    amrex::MultiFab Dx_cumsum( Dx_mf.boxArray(), Dx_mf.DistributionMap(), Dx_mf.nComp(), Dx_mf.nGrow(),
                               amrex::MFInfo() , amrex::FArrayBoxFactory() );
    amrex::MultiFab Dy_cumsum( Dy_mf.boxArray(), Dy_mf.DistributionMap(), Dy_mf.nComp(), Dy_mf.nGrow(),
                               amrex::MFInfo() , amrex::FArrayBoxFactory() );
    amrex::MultiFab Dz_cumsum( Dz_mf.boxArray(), Dz_mf.DistributionMap(), Dz_mf.nComp(), Dz_mf.nGrow(),
                               amrex::MFInfo() , amrex::FArrayBoxFactory() );

    // Compute cumulative sums of D

    // Loop over boxes for Jx
    for (amrex::MFIter mfi(Dx_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dx_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real const> const& Dx_arr = Dx_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dx_cumsum_arr = Dx_cumsum.array( mfi );

        const amrex::Dim3 lo_jx = lbound( Dx_arr );
        const amrex::Dim3 hi_jx = ubound( Dx_arr );

        // Loop over indices within one box
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            Dx_cumsum_arr(i,j,k) = 0.0_rt;
            // Compute cumulative sums of Dx along x
            for ( int ii = lo_jx.x; ii <= i; ++ii ) Dx_cumsum_arr(i,j,k) += Dx_arr(ii,j,k);
        } );
    }

#if (AMREX_SPACEDIM==2) // z direction is along second index

    // Loop over boxes for Jz
    for (amrex::MFIter mfi(Dz_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dz_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real const> const& Dz_arr = Dz_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dz_cumsum_arr = Dz_cumsum.array( mfi );

        const amrex::Dim3 lo_jz = lbound( Dz_arr );
        const amrex::Dim3 hi_jz = ubound( Dz_arr );

        // Loop over indices within one box
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            Dz_cumsum_arr(i,j,k) = 0.0_rt;
            // Compute cumulative sums of Dz along z (second index)
            for ( int jj = lo_jz.y; jj <= j; ++jj ) Dz_cumsum_arr(i,j,k) += Dz_arr(i,jj,k);
        } );
    }

#elif (AMREX_SPACEDIM==3) // z direction is along third index

    // Loop over boxes for Jy
    for (amrex::MFIter mfi(Dy_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dy_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real const> const& Dy_arr = Dy_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dy_cumsum_arr = Dy_cumsum.array( mfi );

        const amrex::Dim3 lo_jy = lbound( Dy_arr );
        const amrex::Dim3 hi_jy = ubound( Dy_arr );

        // Loop over indices within one box
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            Dy_cumsum_arr(i,j,k) = 0.0_rt;
            // Compute cumulative sums of Dy along y
            for ( int jj = lo_jy.y; jj <= j; ++jj ) Dy_cumsum_arr(i,j,k) += Dy_arr(i,jj,k);
        } );
    }

    // Loop over boxes for Jz
    for (amrex::MFIter mfi(Dz_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dz_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real const> const& Dz_arr = Dz_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dz_cumsum_arr = Dz_cumsum.array( mfi );

        const amrex::Dim3 lo_jz = lbound( Dz_arr );
        const amrex::Dim3 hi_jz = ubound( Dz_arr );

        // Loop over indices within one box
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            Dz_cumsum_arr(i,j,k) = 0.0_rt;
            // Compute cumulative sums of Dz along z (third index)
            for ( int kk = lo_jz.z; kk <= k; ++kk ) Dz_cumsum_arr(i,j,k) += Dz_arr(i,j,kk);
        } );
    }

#endif

    // Backward Fourier transform of J
    // (units divided by [L]: one factor [L] from before remains)
    field_data.BackwardTransform( *current[0], Idx::Jx, 0 );
    field_data.BackwardTransform( *current[1], Idx::Jy, 0 );
    field_data.BackwardTransform( *current[2], Idx::Jz, 0 );

    // Subtract from J averages of cumulative sums of D

    // Loop over boxes for Jx
    for (amrex::MFIter mfi(Dx_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dx_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real> const& Dx_arr = Dx_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dx_cumsum_arr = Dx_cumsum.array( mfi );

        const amrex::Dim3 lo_jx = lbound( Dx_arr );
        const amrex::Dim3 hi_jx = ubound( Dx_arr );
        const int nx = hi_jx.x-lo_jx.x+1;

        // Compute average of cumulative sums of Dx along x
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            amrex::Real Dx_cumavg = 0.0_rt;
            for ( int ii = lo_jx.x; ii <= hi_jx.x; ++ii ) Dx_cumavg += Dx_cumsum_arr(ii,j,k);
            Dx_cumavg /= static_cast<Real>(nx);
            Dx_arr(i,j,k) /= dx;
            Dx_arr(i,j,k) -= Dx_cumavg;
        } );
    }

#if (AMREX_SPACEDIM==2) // z direction is along second index

    // Loop over boxes for Jz
    for (amrex::MFIter mfi(Dz_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dz_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real> const& Dz_arr = Dz_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dz_cumsum_arr = Dz_cumsum.array( mfi );

        const amrex::Dim3 lo_jz = lbound( Dz_arr );
        const amrex::Dim3 hi_jz = ubound( Dz_arr );
        const int nz = hi_jz.y-lo_jz.y+1;

        // Compute average of cumulative sums of Dz along z (second index)
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            amrex::Real Dz_cumavg = 0.0_rt;
            for ( int jj = lo_jz.y; jj <= hi_jz.y; ++jj ) Dz_cumavg += Dz_cumsum_arr(i,jj,k);
            Dz_cumavg /= static_cast<Real>(nz);
            Dz_arr(i,j,k) /= dz;
            Dz_arr(i,j,k) -= Dz_cumavg;
        } );
    }

#elif (AMREX_SPACEDIM==3) // z direction is along third index

    // Loop over boxes for Jy
    for (amrex::MFIter mfi(Dy_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dy_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real> const& Dy_arr = Dy_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dy_cumsum_arr = Dy_cumsum.array( mfi );

        const amrex::Dim3 lo_jy = lbound( Dy_arr );
        const amrex::Dim3 hi_jy = ubound( Dy_arr );
        const int ny = hi_jy.y-lo_jy.y+1;

        // Compute average of cumulative sums of Dy along y
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            amrex::Real Dy_cumavg = 0.0_rt;
            for ( int jj = lo_jy.y; jj <= hi_jy.y; ++jj ) Dy_cumavg += Dy_cumsum_arr(i,jj,k);
            Dy_cumavg /= static_cast<Real>(ny);
            Dy_arr(i,j,k) /= dy;
            Dy_arr(i,j,k) -= Dy_cumavg;
        } );
    }

    // Loop over boxes for Jz
    for (amrex::MFIter mfi(Dz_mf); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = Dz_mf[mfi].box();

        // Original current D deposited in CurrentDeposition.H
        amrex::Array4<amrex::Real> const& Dz_arr = Dz_mf.array( mfi );

        // 3D array to store directional cumulative sum of D
        amrex::Array4<amrex::Real> const Dz_cumsum_arr = Dz_cumsum.array( mfi );

        const amrex::Dim3 lo_jz = lbound( Dz_arr );
        const amrex::Dim3 hi_jz = ubound( Dz_arr );
        const int nz = hi_jz.z-lo_jz.z+1;

        // Compute average of cumulative sums of Dz along z (third index)
        ParallelFor( bx, [=] AMREX_GPU_DEVICE( int i, int j, int k ) noexcept
        {
            amrex::Real Dz_cumavg = 0.0_rt;
            for ( int kk = lo_jz.z; kk <= hi_jz.z; ++kk ) Dz_cumavg += Dz_cumsum_arr(i,j,kk);
            Dz_cumavg /= static_cast<Real>(nz);
            Dz_arr(i,j,k) /= dz;
            Dz_arr(i,j,k) -= Dz_cumavg;
        } );
    }

#endif
}
#endif // WARPX_USE_PSATD
