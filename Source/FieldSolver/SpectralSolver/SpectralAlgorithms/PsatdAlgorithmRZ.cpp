#include <PsatdAlgorithmRZ.H>
#include <WarpXConst.H>
#include <cmath>

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
PsatdAlgorithmRZ::PsatdAlgorithmRZ (SpectralHankelKSpace const & spectral_kspace,
                                    DistributionMapping const & dm,
                                    int const n_rz_azimuthal_modes, int const norder_z,
                                    bool const nodal, Real const dt)
     // Initialize members of base class
     : SpectralBaseAlgorithmRZ(spectral_kspace, dm,
                               norder_z, nodal),
       dt(dt)
{

    // Allocate the arrays of coefficients
    BoxArray const & ba = spectral_kspace.spectralspace_ba;
    C_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    S_ck_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    X1_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    X2_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    X3_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);

    coefficients_initialized = false;
}

/* Advance the E and B field in spectral space (stored in `f`)
 * over one time step */
void
PsatdAlgorithmRZ::pushSpectralFields(SpectralFieldDataHankel & f)
{

    if (not coefficients_initialized) {
        // This is called from here since it needs the kr values
        // which can be obtained from the SpectralFieldDataHankel
        InitializeSpectralCoefficients(f);
        coefficients_initialized = true;
    }

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        Box const & bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();
        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
        Array4<const Real> X1_arr = X1_coef[mfi].array();
        Array4<const Real> X2_arr = X2_coef[mfi].array();
        Array4<const Real> X3_arr = X3_coef[mfi].array();

        // Extract pointers for the k vectors
        auto const & kr = f.getKrArray(mfi);
        Real const* kr_arr = kr.dataPtr();
        Real const* modified_kz_arr = modified_kz_vec[mfi].dataPtr();
        int const nr = bx.length(0);

        // Loop over indices within one box
        // Note that k = 0
        int const modes = f.n_rz_azimuthal_modes;
        ParallelFor(bx, modes,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept
        {

            using Idx = SpectralFieldIndex;
            auto const Ep_m = Idx::Ex + Idx::n_fields*mode;
            auto const Em_m = Idx::Ey + Idx::n_fields*mode;
            auto const Ez_m = Idx::Ez + Idx::n_fields*mode;
            auto const Bp_m = Idx::Bx + Idx::n_fields*mode;
            auto const Bm_m = Idx::By + Idx::n_fields*mode;
            auto const Bz_m = Idx::Bz + Idx::n_fields*mode;

            // Record old values of the fields to be updated
            Complex const Ep_old = fields(i,j,k,Ep_m);
            Complex const Em_old = fields(i,j,k,Em_m);
            Complex const Ez_old = fields(i,j,k,Ez_m);
            Complex const Bp_old = fields(i,j,k,Bp_m);
            Complex const Bm_old = fields(i,j,k,Bm_m);
            Complex const Bz_old = fields(i,j,k,Bz_m);
            // Shortcut for the values of J and rho
            Complex const Jp = fields(i,j,k,Idx::Jx + Idx::n_fields*mode);
            Complex const Jm = fields(i,j,k,Idx::Jy + Idx::n_fields*mode);
            Complex const Jz = fields(i,j,k,Idx::Jz + Idx::n_fields*mode);
            Complex const rho_old = fields(i,j,k,Idx::rho_old + Idx::n_fields*mode);
            Complex const rho_new = fields(i,j,k,Idx::rho_new + Idx::n_fields*mode);

            // k vector values, and coefficients
            int const ii = i + nr*mode;
            Real const kr = kr_arr[ii];
            Real const kz = modified_kz_arr[j];

            constexpr Real c2 = PhysConst::c*PhysConst::c;
            constexpr Real inv_ep0 = 1./PhysConst::ep0;
            Complex const I = Complex{0,1};
            Real const C = C_arr(i,j,k,mode);
            Real const S_ck = S_ck_arr(i,j,k,mode);
            Real const X1 = X1_arr(i,j,k,mode);
            Real const X2 = X2_arr(i,j,k,mode);
            Real const X3 = X3_arr(i,j,k,mode);

            // Update E (see WarpX online documentation: theory section)
            fields(i,j,k,Ep_m) = C*Ep_old
                        + S_ck*(-c2*I*kr/2.*Bz_old + c2*kz*Bp_old - inv_ep0*Jp)
                        + kr*(X2*rho_new - X3*rho_old);
            fields(i,j,k,Em_m) = C*Em_old
                        + S_ck*(-c2*I*kr/2.*Bz_old - c2*kz*Bm_old - inv_ep0*Jm)
                        - kr*(X2*rho_new - X3*rho_old);
            fields(i,j,k,Ez_m) = C*Ez_old
                        + S_ck*(c2*I*kr*Bp_old + c2*I*kr*Bm_old - inv_ep0*Jz)
                        - I*kz*(X2*rho_new - X3*rho_old);
            // Update B (see WarpX online documentation: theory section)
            fields(i,j,k,Bp_m) = C*Bp_old
                        - S_ck*(-I*kr/2.*Ez_old + kz*Ep_old)
                        + X1*(-I*kr/2.*Jz + kz*Jp);
            fields(i,j,k,Bm_m) = C*Bm_old
                        - S_ck*(-I*kr/2.*Ez_old - kz*Em_old)
                        + X1*(-I*kr/2.*Jz - kz*Jm);
            fields(i,j,k,Bz_m) = C*Bz_old
                        - S_ck*I*(kr*Ep_old + kr*Em_old)
                        + X1*I*(kr*Jp + kr*Jm);
        });
    }
};

void PsatdAlgorithmRZ::InitializeSpectralCoefficients (SpectralFieldDataHankel const & f)
{

    // Fill them with the right values:
    // Loop over boxes and allocate the corresponding coefficients
    // for each box owned by the local MPI proc
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        Box const & bx = f.fields[mfi].box();

        // Extract pointers for the k vectors
        Real* const modified_kz = modified_kz_vec[mfi].dataPtr();

        // Extract arrays for the coefficients
        Array4<Real> C = C_coef[mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();
        Array4<Real> X1 = X1_coef[mfi].array();
        Array4<Real> X2 = X2_coef[mfi].array();
        Array4<Real> X3 = X3_coef[mfi].array();

        auto const & kr = f.getKrArray(mfi);
        Real const* kr_arr = kr.dataPtr();
        int const nr = bx.length(0);

        // Loop over indices within one box
        int const modes = f.n_rz_azimuthal_modes;
        ParallelFor(bx, modes,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept
        {
            // Calculate norm of vector
            int const ii = i + nr*mode;
            Real const kr = kr_arr[ii];
            Real const kz = modified_kz[j];
            Real const k_norm = std::sqrt(kr*kr + kz*kz);

            // Calculate coefficients
            constexpr Real c = PhysConst::c;
            constexpr Real ep0 = PhysConst::ep0;
            if (k_norm != 0){
                C(i,j,k,mode) = std::cos(c*k_norm*dt);
                S_ck(i,j,k,mode) = std::sin(c*k_norm*dt)/(c*k_norm);
                X1(i,j,k,mode) = (1. - C(i,j,k,mode))/(ep0 * c*c * k_norm*k_norm);
                X2(i,j,k,mode) = (1. - S_ck(i,j,k,mode)/dt)/(ep0 * k_norm*k_norm);
                X3(i,j,k,mode) = (C(i,j,k,mode) - S_ck(i,j,k,mode)/dt)/(ep0 * k_norm*k_norm);
            } else { // Handle k_norm = 0, by using the analytical limit
                C(i,j,k,mode) = 1.;
                S_ck(i,j,k,mode) = dt;
                X1(i,j,k,mode) = 0.5 * dt*dt / ep0;
                X2(i,j,k,mode) = c*c * dt*dt / (6.*ep0);
                X3(i,j,k,mode) = - c*c * dt*dt / (3.*ep0);
            }
        });
     }
}
