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
                         const int norder_z, const bool nodal, const Real dt,
                         const bool update_with_rho)
    // Initialize members of base class
    : m_dt( dt ),
      m_update_with_rho( update_with_rho ),
      SpectralBaseAlgorithm( spectral_kspace, dm, norder_x, norder_y, norder_z, nodal )
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
}

/**
 * \brief Advance E and B fields in spectral space (stored in `f`) over one time step
 */
void
PsatdAlgorithm::pushSpectralFields(SpectralFieldData& f) const{

    const bool update_with_rho = m_update_with_rho;

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
        ParallelFor( bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
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
            constexpr Real inv_eps0 = 1.0_rt/PhysConst::ep0;

            const Complex I = Complex{0,1};

            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Real X1 = X1_arr(i,j,k);
            const Real X2 = X2_arr(i,j,k);
            const Real X3 = X3_arr(i,j,k);

            // Update E (see WarpX online documentation: theory section)

            if (update_with_rho) {

                fields(i,j,k,Idx::Ex) = C*Ex_old + S_ck*(c2*I*(ky*Bz_old-kz*By_old)-inv_eps0*Jx)
                                        - I*(X2*rho_new-X3*rho_old)*kx;

                fields(i,j,k,Idx::Ey) = C*Ey_old + S_ck*(c2*I*(kz*Bx_old-kx*Bz_old)-inv_eps0*Jy)
                                        - I*(X2*rho_new-X3*rho_old)*ky;

                fields(i,j,k,Idx::Ez) = C*Ez_old + S_ck*(c2*I*(kx*By_old-ky*Bx_old)-inv_eps0*Jz)
                                        - I*(X2*rho_new-X3*rho_old)*kz;
            } else {

                Complex k_dot_J = kx*Jx + ky*Jy + kz*Jz;
                Complex k_dot_E = kx*Ex_old + ky*Ey_old + kz*Ez_old;

                fields(i,j,k,Idx::Ex) = C*Ex_old + S_ck*(c2*I*(ky*Bz_old-kz*By_old)-inv_eps0*Jx)
                                        + X2*k_dot_E*kx + X3*inv_eps0*k_dot_J*kx;

                fields(i,j,k,Idx::Ey) = C*Ey_old + S_ck*(c2*I*(kz*Bx_old-kx*Bz_old)-inv_eps0*Jy)
                                        + X2*k_dot_E*ky + X3*inv_eps0*k_dot_J*ky;

                fields(i,j,k,Idx::Ez) = C*Ez_old + S_ck*(c2*I*(kx*By_old-ky*Bx_old)-inv_eps0*Jz)
                                        + X2*k_dot_E*kz + X3*inv_eps0*k_dot_J*kz;
            }

            // Update B (see WarpX online documentation: theory section)

            fields(i,j,k,Idx::Bx) = C*Bx_old - S_ck*I*(ky*Ez_old-kz*Ey_old) + X1*I*(ky*Jz-kz*Jy);

            fields(i,j,k,Idx::By) = C*By_old - S_ck*I*(kz*Ex_old-kx*Ez_old) + X1*I*(kz*Jx-kx*Jz);

            fields(i,j,k,Idx::Bz) = C*Bz_old - S_ck*I*(kx*Ey_old-ky*Ex_old) + X1*I*(kx*Jy-ky*Jx);
        } );
    }
};

/**
 * \brief Initialize coefficients for update equations
 */
void PsatdAlgorithm::InitializeSpectralCoefficients(const SpectralKSpace& spectral_kspace,
                                    const amrex::DistributionMapping& dm,
                                    const amrex::Real dt)
{
    const bool update_with_rho = m_update_with_rho;

    const BoxArray& ba = spectral_kspace.spectralspace_ba;

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
        ParallelFor( bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of vector
            const Real k_norm = std::sqrt(
                std::pow(modified_kx[i],2) +
#if (AMREX_SPACEDIM==3)
                std::pow(modified_ky[j],2) + std::pow(modified_kz[k],2));
#else
                std::pow(modified_kz[j],2));
#endif
            // Calculate coefficients
            constexpr Real c = PhysConst::c;
            constexpr Real eps0 = PhysConst::ep0;

            if (k_norm != 0) {
                C(i,j,k)    = std::cos(c*k_norm*dt);
                S_ck(i,j,k) = std::sin(c*k_norm*dt)/(c*k_norm);
                X1(i,j,k) = (1.0_rt-C(i,j,k))/(eps0*c*c*k_norm*k_norm);
                if (update_with_rho) {
                    X2(i,j,k) = (1.0_rt-S_ck(i,j,k)/dt)/(eps0*k_norm*k_norm);
                    X3(i,j,k) = (C(i,j,k)-S_ck(i,j,k)/dt)/(eps0*k_norm*k_norm);
                } else {
                    X2(i,j,k) = (1.0_rt-C(i,j,k))/(k_norm*k_norm);
                    X3(i,j,k) = (S_ck(i,j,k)-dt)/(k_norm*k_norm);
                }
            } else { // Handle k_norm = 0 with analytical limit
                C(i,j,k) = 1.0_rt;
                S_ck(i,j,k) = dt;
                X1(i,j,k) = 0.5_rt*dt*dt/eps0;
                if (update_with_rho) {
                    X2(i,j,k) = c*c*dt*dt/(6.0_rt*eps0);
                    X3(i,j,k) = -c*c*dt*dt/(3.0_rt*eps0);
                } else {
                    X2(i,j,k) = 0.5_rt*dt*dt*c*c;
                    X3(i,j,k) = -c*c*dt*dt*dt/6.0_rt;
                }
            }
        } );
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
        ParallelFor( bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
            const Real k_norm = std::sqrt( kx*kx + ky*ky + kz*kz );

            constexpr Complex I = Complex{0,1};

            // div(J) in Fourier space
            const Complex k_dot_J = kx*Jx + ky*Jy + kz*Jz;

            // Correct J
            if ( k_norm != 0 )
            {
                fields(i,j,k,Idx::Jx) = Jx - (k_dot_J-I*(rho_new-rho_old)/dt)*kx/(k_norm*k_norm);
                fields(i,j,k,Idx::Jy) = Jy - (k_dot_J-I*(rho_new-rho_old)/dt)*ky/(k_norm*k_norm);
                fields(i,j,k,Idx::Jz) = Jz - (k_dot_J-I*(rho_new-rho_old)/dt)*kz/(k_norm*k_norm);
            }
        } );
    }

    // Backward Fourier transform of J
    field_data.BackwardTransform( *current[0], Idx::Jx, 0 );
    field_data.BackwardTransform( *current[1], Idx::Jy, 0 );
    field_data.BackwardTransform( *current[2], Idx::Jz, 0 );
}
#endif // WARPX_USE_PSATD
