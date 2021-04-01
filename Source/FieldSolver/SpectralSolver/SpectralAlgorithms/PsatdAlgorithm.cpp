/* Copyright 2019
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PsatdAlgorithm.H"
#include "Utils/WarpXConst.H"

#include <cmath>

#if WARPX_USE_PSATD

using namespace amrex;

PsatdAlgorithm::PsatdAlgorithm(
    const SpectralKSpace& spectral_kspace,
    const DistributionMapping& dm,
    const int norder_x,
    const int norder_y,
    const int norder_z,
    const bool nodal,
    const amrex::Array<amrex::Real,3>& v_galilean,
    const amrex::Real dt,
    const bool update_with_rho,
    const bool time_averaging)
    // Initializer list
    : SpectralBaseAlgorithm(spectral_kspace, dm, norder_x, norder_y, norder_z, nodal),
    // Initialize the centered finite-order modified k vectors:
    // these are computed always with the assumption of centered grids
    // (argument nodal = true), for both nodal and staggered simulations
    modified_kx_vec_centered(spectral_kspace.getModifiedKComponent(dm, 0, norder_x, true)),
#if (AMREX_SPACEDIM == 3)
    modified_ky_vec_centered(spectral_kspace.getModifiedKComponent(dm, 1, norder_y, true)),
    modified_kz_vec_centered(spectral_kspace.getModifiedKComponent(dm, 2, norder_z, true)),
#else
    modified_kz_vec_centered(spectral_kspace.getModifiedKComponent(dm, 1, norder_z, true)),
#endif
    m_v_galilean(v_galilean),
    m_dt(dt),
    m_update_with_rho(update_with_rho),
    m_time_averaging(time_averaging)
{
    const amrex::BoxArray& ba = spectral_kspace.spectralspace_ba;

    m_is_galilean = (v_galilean[0] != 0.) || (v_galilean[1] != 0.) || (v_galilean[2] != 0.);

    // Always allocate these coefficients
    C_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X3_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

    if (m_is_galilean)
    {
        X4_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        T2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    }

    InitializeSpectralCoefficients(spectral_kspace, dm, dt);

    // Allocate these coefficients only with averaged Galilean PSATD
    if (time_averaging)
    {
        Psi1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Psi2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Y1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Y3_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Y2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Y4_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

        InitializeSpectralCoefficientsAveraging(spectral_kspace, dm, dt);
    }
}

void
PsatdAlgorithm::pushSpectralFields (SpectralFieldData& f) const
{
    const bool update_with_rho = m_update_with_rho;
    const bool time_averaging  = m_time_averaging;
    const bool is_galilean     = m_is_galilean;

    // Loop over boxes
    for (amrex::MFIter mfi(f.fields); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> fields = f.fields[mfi].array();

        // These coefficients are always allocated
        amrex::Array4<const amrex::Real> C_arr = C_coef[mfi].array();
        amrex::Array4<const amrex::Real> S_ck_arr = S_ck_coef[mfi].array();
        amrex::Array4<const Complex> X1_arr = X1_coef[mfi].array();
        amrex::Array4<const Complex> X2_arr = X2_coef[mfi].array();
        amrex::Array4<const Complex> X3_arr = X3_coef[mfi].array();

        amrex::Array4<const Complex> X4_arr;
        amrex::Array4<const Complex> T2_arr;
        if (is_galilean)
        {
            X4_arr = X4_coef[mfi].array();
            T2_arr = T2_coef[mfi].array();
        }

        // These coefficients are allocated only with averaged Galilean PSATD
        amrex::Array4<const Complex> Psi1_arr;
        amrex::Array4<const Complex> Psi2_arr;
        amrex::Array4<const Complex> Y1_arr;
        amrex::Array4<const Complex> Y2_arr;
        amrex::Array4<const Complex> Y3_arr;
        amrex::Array4<const Complex> Y4_arr;

        if (time_averaging)
        {
            Psi1_arr = Psi1_coef[mfi].array();
            Psi2_arr = Psi2_coef[mfi].array();
            Y1_arr = Y1_coef[mfi].array();
            Y2_arr = Y2_coef[mfi].array();
            Y3_arr = Y3_coef[mfi].array();
            Y4_arr = Y4_coef[mfi].array();
        }

        // Extract pointers for the k vectors
        const amrex::Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
        const amrex::Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            using Idx = SpectralFieldIndex;
            using AvgIdx = SpectralAvgFieldIndex;

            // Record old values of the fields to be updated
            const Complex Ex_old = fields(i,j,k,Idx::Ex);
            const Complex Ey_old = fields(i,j,k,Idx::Ey);
            const Complex Ez_old = fields(i,j,k,Idx::Ez);
            const Complex Bx_old = fields(i,j,k,Idx::Bx);
            const Complex By_old = fields(i,j,k,Idx::By);
            const Complex Bz_old = fields(i,j,k,Idx::Bz);

            // Shortcuts for the values of J and rho
            const Complex Jx = fields(i,j,k,Idx::Jx);
            const Complex Jy = fields(i,j,k,Idx::Jy);
            const Complex Jz = fields(i,j,k,Idx::Jz);
            const Complex rho_old = fields(i,j,k,Idx::rho_old);
            const Complex rho_new = fields(i,j,k,Idx::rho_new);

            // k vector values
            const amrex::Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM == 3)
            const amrex::Real ky = modified_ky_arr[j];
            const amrex::Real kz = modified_kz_arr[k];
#else
            constexpr amrex::Real ky = 0._rt;
            const     amrex::Real kz = modified_kz_arr[j];
#endif
            // Physical constants and imaginary unit
            const amrex::Real c2 = std::pow(PhysConst::c, 2);
            constexpr Complex I = Complex{0._rt, 1._rt};

            // These coefficients are initialized in the function InitializeSpectralCoefficients
            const amrex::Real C = C_arr(i,j,k);
            const amrex::Real S_ck = S_ck_arr(i,j,k);
            const Complex X1 = X1_arr(i,j,k);
            const Complex X2 = X2_arr(i,j,k);
            const Complex X3 = X3_arr(i,j,k);
            const Complex X4 = (is_galilean) ? X4_arr(i,j,k) : - S_ck / PhysConst::ep0;
            const Complex T2 = (is_galilean) ? T2_arr(i,j,k) : 1.0_rt;

            // Update equations for E in the formulation with rho
            // T2 = 1 always with standard PSATD (zero Galilean velocity)

            if (update_with_rho)
            {
                fields(i,j,k,Idx::Ex) = T2 * C * Ex_old
                                        + I * c2 * T2 * S_ck * (ky * Bz_old - kz * By_old)
                                        + X4 * Jx - I * (X2 * rho_new - T2 * X3 * rho_old) * kx;

                fields(i,j,k,Idx::Ey) = T2 * C * Ey_old
                                        + I * c2 * T2 * S_ck * (kz * Bx_old - kx * Bz_old)
                                        + X4 * Jy - I * (X2 * rho_new - T2 * X3 * rho_old) * ky;

                fields(i,j,k,Idx::Ez) = T2 * C * Ez_old
                                        + I * c2 * T2 * S_ck * (kx * By_old - ky * Bx_old)
                                        + X4 * Jz - I * (X2 * rho_new - T2 * X3 * rho_old) * kz;
            }

            // Update equations for E in the formulation without rho
            // T2 = 1 always with standard PSATD (zero Galilean velocity)

            else {

                Complex k_dot_J = kx * Jx + ky * Jy + kz * Jz;
                Complex k_dot_E = kx * Ex_old + ky * Ey_old + kz * Ez_old;

                fields(i,j,k,Idx::Ex) = T2 * C * Ex_old
                                        + I * c2 * T2 * S_ck * (ky * Bz_old - kz * By_old)
                                        + X4 * Jx + X2 * k_dot_E * kx + X3 * k_dot_J * kx;

                fields(i,j,k,Idx::Ey) = T2 * C * Ey_old
                                        + I * c2 * T2 * S_ck * (kz * Bx_old - kx * Bz_old)
                                        + X4 * Jy + X2 * k_dot_E * ky + X3 * k_dot_J * ky;

                fields(i,j,k,Idx::Ez) = T2 * C * Ez_old
                                        + I * c2 * T2 * S_ck * (kx * By_old - ky * Bx_old)
                                        + X4 * Jz + X2 * k_dot_E * kz + X3 * k_dot_J * kz;
            }

            // Update equations for B
            // T2 = 1 always with standard PSATD (zero Galilean velocity)

            fields(i,j,k,Idx::Bx) = T2 * C * Bx_old
                                    - I * T2 * S_ck * (ky * Ez_old - kz * Ey_old)
                                    + I * X1 * (ky * Jz - kz * Jy);

            fields(i,j,k,Idx::By) = T2 * C * By_old
                                    - I * T2 * S_ck * (kz * Ex_old - kx * Ez_old)
                                    + I * X1 * (kz * Jx - kx * Jz);

            fields(i,j,k,Idx::Bz) = T2 * C * Bz_old
                                    - I * T2 * S_ck * (kx * Ey_old - ky * Ex_old)
                                    + I * X1 * (kx * Jy - ky * Jx);

            // Additional update equations for averaged Galilean algorithm

            if (time_averaging)
            {
                // These coefficients are initialized in the function InitializeSpectralCoefficients below
                const Complex Psi1 = Psi1_arr(i,j,k);
                const Complex Psi2 = Psi2_arr(i,j,k);
                const Complex Y1 = Y1_arr(i,j,k);
                const Complex Y3 = Y3_arr(i,j,k);
                const Complex Y2 = Y2_arr(i,j,k);
                const Complex Y4 = Y4_arr(i,j,k);

                fields(i,j,k,AvgIdx::Ex_avg) = Psi1 * Ex_old
                                               - I * c2 * Psi2 * (ky * Bz_old - kz * By_old)
                                               + Y4 * Jx + (Y2 * rho_new + Y3 * rho_old) * kx;

                fields(i,j,k,AvgIdx::Ey_avg) = Psi1 * Ey_old
                                               - I * c2 * Psi2 * (kz * Bx_old - kx * Bz_old)
                                               + Y4 * Jy + (Y2 * rho_new + Y3 * rho_old) * ky;

                fields(i,j,k,AvgIdx::Ez_avg) = Psi1 * Ez_old
                                               - I * c2 * Psi2 * (kx * By_old - ky * Bx_old)
                                               + Y4 * Jz + (Y2 * rho_new + Y3 * rho_old) * kz;

                fields(i,j,k,AvgIdx::Bx_avg) = Psi1 * Bx_old
                                               + I * Psi2 * (ky * Ez_old - kz * Ey_old)
                                               + I * Y1 * (ky * Jz - kz * Jy);

                fields(i,j,k,AvgIdx::By_avg) = Psi1 * By_old
                                               + I * Psi2 * (kz * Ex_old - kx * Ez_old)
                                               + I * Y1 * (kz * Jx - kx * Jz);

                fields(i,j,k,AvgIdx::Bz_avg) = Psi1 * Bz_old
                                               + I * Psi2 * (kx * Ey_old - ky * Ex_old)
                                               + I * Y1 * (kx * Jy - ky * Jx);
            }
        });
    }
}

void PsatdAlgorithm::InitializeSpectralCoefficients (
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const amrex::Real dt)
{
    const bool update_with_rho = m_update_with_rho;
    const bool is_galilean     = m_is_galilean;

    const amrex::BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Loop over boxes and allocate the corresponding coefficients for each box
    for (amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const amrex::Real* kx_s = modified_kx_vec[mfi].dataPtr();
        const amrex::Real* kx_c = modified_kx_vec_centered[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
        const amrex::Real* ky_s = modified_ky_vec[mfi].dataPtr();
        const amrex::Real* ky_c = modified_ky_vec_centered[mfi].dataPtr();
#endif
        const amrex::Real* kz_s = modified_kz_vec[mfi].dataPtr();
        const amrex::Real* kz_c = modified_kz_vec_centered[mfi].dataPtr();

        // Coefficients always allocated
        amrex::Array4<amrex::Real> C = C_coef[mfi].array();
        amrex::Array4<amrex::Real> S_ck = S_ck_coef[mfi].array();
        amrex::Array4<Complex> X1 = X1_coef[mfi].array();
        amrex::Array4<Complex> X2 = X2_coef[mfi].array();
        amrex::Array4<Complex> X3 = X3_coef[mfi].array();

        amrex::Array4<Complex> X4;
        amrex::Array4<Complex> T2;
        if (is_galilean)
        {
            X4 = X4_coef[mfi].array();
            T2 = T2_coef[mfi].array();
        }

        // Extract Galilean velocity
        amrex::Real vg_x = m_v_galilean[0];
#if (AMREX_SPACEDIM == 3)
        amrex::Real vg_y = m_v_galilean[1];
#endif
        amrex::Real vg_z = m_v_galilean[2];

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of k vector
            const amrex::Real knorm_s = std::sqrt(
                std::pow(kx_s[i], 2) +
#if (AMREX_SPACEDIM == 3)
                std::pow(ky_s[j], 2) + std::pow(kz_s[k], 2));
#else
                std::pow(kz_s[j], 2));
#endif
            // Physical constants and imaginary unit
            constexpr amrex::Real c = PhysConst::c;
            constexpr amrex::Real ep0 = PhysConst::ep0;
            constexpr Complex I = Complex{0._rt, 1._rt};

            const amrex::Real c2 = std::pow(c, 2);
            const amrex::Real dt2 = std::pow(dt, 2);
            const amrex::Real dt3 = std::pow(dt, 3);

            // Calculate the dot product of the k vector with the Galilean velocity.
            // This has to be computed always with the centered (that is, nodal) finite-order
            // modified k vectors, to work correctly for both nodal and staggered simulations.
            // w_c = 0 always with standard PSATD (zero Galilean velocity).
            const amrex::Real w_c = kx_c[i]*vg_x +
#if (AMREX_SPACEDIM == 3)
                ky_c[j]*vg_y + kz_c[k]*vg_z;
#else
                kz_c[j]*vg_z;
#endif
            const amrex::Real w2_c = std::pow(w_c, 2);

            const amrex::Real om_s = c * knorm_s;
            const amrex::Real om2_s = std::pow(om_s, 2);

            const Complex theta_c      = amrex::exp( I * w_c * dt * 0.5_rt);
            const Complex theta2_c     = amrex::exp( I * w_c * dt);
            const Complex theta_c_star = amrex::exp(-I * w_c * dt * 0.5_rt);

            // C
            C(i,j,k) = std::cos(om_s * dt);

            // S_ck
            if (om_s != 0.)
            {
                S_ck(i,j,k) = std::sin(om_s * dt) / om_s;
            }
            else // om_s = 0
            {
                S_ck(i,j,k) = dt;
            }

            // Auxiliary variable
            amrex::Real tmp;
            if (om_s != 0.)
            {
                tmp = (1._rt - C(i,j,k)) / (ep0 * om2_s);
            }
            else // om_s = 0
            {
                tmp = 0.5_rt * dt2 / ep0;
            }

            // T2
            if (is_galilean)
            {
                T2(i,j,k) = theta_c * theta_c;
            }

            // X1 (multiplies i*([k] \times J) in the update equation for update B)
            if ((om_s != 0.) || (w_c != 0.))
            {
                X1(i,j,k) = (1._rt - theta2_c * C(i,j,k) + I * w_c * theta2_c * S_ck(i,j,k))
                            / (ep0 * (om2_s - w2_c));
            }
            else // om_s = 0 and w_c = 0
            {
                X1(i,j,k) = 0.5_rt * dt2 / ep0;
            }

            // X2 (multiplies rho_new      if update_with_rho = 1 in the update equation for E)
            // X2 (multiplies ([k] \dot E) if update_with_rho = 0 in the update equation for E)
            if (update_with_rho)
            {
                if (w_c != 0.)
                {
                    X2(i,j,k) = c2 * (theta_c_star * X1(i,j,k) - theta_c * tmp)
                                / (theta_c_star - theta_c);
                }
                else // w_c = 0
                {
                    if (om_s != 0.)
                    {
                        X2(i,j,k) = c2 * (dt - S_ck(i,j,k)) / (ep0 * dt * om2_s);
                    }
                    else // om_s = 0 and w_c = 0
                    {
                        X2(i,j,k) = c2 * dt2 / (6._rt * ep0);
                    }
                }
            }
            else // update_with_rho = 0
            {
                X2(i,j,k) = c2 * ep0 * theta2_c * tmp;
            }

            // X3 (multiplies rho_old      if update_with_rho = 1 in the update equation for E)
            // X3 (multiplies ([k] \dot J) if update_with_rho = 0 in the update equation for E)
            if (update_with_rho)
            {
                if (w_c != 0.)
                {
                    X3(i,j,k) = c2 * (theta_c_star * X1(i,j,k) - theta_c_star * tmp)
                                / (theta_c_star - theta_c);
                }
                else // w_c = 0
                {
                    if (om_s != 0.)
                    {
                        X3(i,j,k) = c2 * (dt * C(i,j,k) - S_ck(i,j,k)) / (ep0 * dt * om2_s);
                    }
                    else // om_s = 0 and w_c = 0
                    {
                        X3(i,j,k) = - c2 * dt2 / (3._rt * ep0);
                    }
                }
            }
            else // update_with_rho = 0
            {
                if (w_c != 0.)
                {
                    X3(i,j,k) = I * c2 * (theta2_c * tmp - X1(i,j,k)) / w_c;
                }
                else // w_c = 0
                {
                    if (om_s != 0.)
                    {
                        X3(i,j,k) = c2 * (S_ck(i,j,k) - dt) / (ep0 * om2_s);
                    }
                    else // om_s = 0 and w_c = 0
                    {
                        X3(i,j,k) = - c2 * dt3 / (6._rt * ep0);
                    }
                }
            }

            // X4 (multiplies J in the update equation for E)
            if (is_galilean)
            {
                X4(i,j,k) = I * w_c * X1(i,j,k) - theta2_c * S_ck(i,j,k) / ep0;
            }
        });
    }
}

void PsatdAlgorithm::InitializeSpectralCoefficientsAveraging (
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const amrex::Real dt)
{
    const amrex::BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Loop over boxes and allocate the corresponding coefficients for each box
    for (amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const amrex::Real* kx_s = modified_kx_vec[mfi].dataPtr();
        const amrex::Real* kx_c = modified_kx_vec_centered[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
        const amrex::Real* ky_s = modified_ky_vec[mfi].dataPtr();
        const amrex::Real* ky_c = modified_ky_vec_centered[mfi].dataPtr();
#endif
        const amrex::Real* kz_s = modified_kz_vec[mfi].dataPtr();
        const amrex::Real* kz_c = modified_kz_vec_centered[mfi].dataPtr();

        // Coefficients allocated only with averaged Galilean PSATD
        amrex::Array4<Complex> Psi1 = Psi1_coef[mfi].array();
        amrex::Array4<Complex> Psi2 = Psi2_coef[mfi].array();
        amrex::Array4<Complex> Y1 = Y1_coef[mfi].array();
        amrex::Array4<Complex> Y3 = Y3_coef[mfi].array();
        amrex::Array4<Complex> Y2 = Y2_coef[mfi].array();
        amrex::Array4<Complex> Y4 = Y4_coef[mfi].array();

        // Extract Galilean velocity
        amrex::Real vg_x = m_v_galilean[0];
#if (AMREX_SPACEDIM == 3)
        amrex::Real vg_y = m_v_galilean[1];
#endif
        amrex::Real vg_z = m_v_galilean[2];

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of k vector
            const amrex::Real knorm_s = std::sqrt(
                std::pow(kx_s[i], 2) +
#if (AMREX_SPACEDIM == 3)
                std::pow(ky_s[j], 2) + std::pow(kz_s[k], 2));
#else
                std::pow(kz_s[j], 2));
#endif
            // Physical constants and imaginary unit
            constexpr amrex::Real c = PhysConst::c;
            constexpr amrex::Real ep0 = PhysConst::ep0;
            constexpr Complex I = Complex{0._rt, 1._rt};

            const amrex::Real c2 = std::pow(c, 2);
            const amrex::Real dt2 = std::pow(dt, 2);

            // Calculate the dot product of the k vector with the Galilean velocity.
            // This has to be computed always with the centered (that is, nodal) finite-order
            // modified k vectors, to work correctly for both nodal and staggered simulations.
            // w_c = 0 always with standard PSATD (zero Galilean velocity).
            const amrex::Real w_c = kx_c[i]*vg_x +
#if (AMREX_SPACEDIM == 3)
                ky_c[j]*vg_y + kz_c[k]*vg_z;
#else
                kz_c[j]*vg_z;
#endif
            const amrex::Real w2_c = std::pow(w_c, 2);
            const amrex::Real w3_c = std::pow(w_c, 3);

            const amrex::Real om_s = c * knorm_s;
            const amrex::Real om2_s = std::pow(om_s, 2);
            const amrex::Real om4_s = std::pow(om_s, 4);

            const Complex theta_c  = amrex::exp(I * w_c * dt * 0.5_rt);
            const Complex theta2_c = amrex::exp(I * w_c * dt);
            const Complex theta3_c = amrex::exp(I * w_c * dt * 1.5_rt);
            const Complex theta5_c = amrex::exp(I * w_c * dt * 2.5_rt);

            // C1,C3
            const amrex::Real C1 = std::cos(0.5_rt * om_s * dt);
            const amrex::Real C3 = std::cos(1.5_rt * om_s * dt);

            // S1_om, S3_om
            amrex::Real S1_om, S3_om;
            if (om_s != 0.)
            {
                S1_om = std::sin(0.5_rt * om_s * dt) / om_s;
                S3_om = std::sin(1.5_rt * om_s * dt) / om_s;
            }
            else // om_s = 0
            {
                S1_om = 0.5_rt * dt;
                S3_om = 1.5_rt * dt;
            }

            // Psi1 (multiplies E in the update equation for <E>)
            // Psi1 (multiplies B in the update equation for <B>)
            if ((om_s != 0.) || (w_c != 0.))
            {
                Psi1(i,j,k) = (theta3_c * (om2_s * S3_om + I * w_c * C3)
                              - theta_c * (om2_s * S1_om + I * w_c * C1)) / (dt * (om2_s - w2_c));
            }
            else // om_s = 0 and w_c = 0
            {
                Psi1(i,j,k) = 1._rt;
            }

            // Psi2 (multiplies i*([k] \times B) in the update equation for <E>)
            // Psi2 (multiplies i*([k] \times E) in the update equation for <B>)
            if ((om_s != 0.) || (w_c != 0.))
            {
                Psi2(i,j,k) = (theta3_c * (C3 - I * w_c * S3_om)
                              - theta_c * (C1 - I * w_c * S1_om)) / (dt * (om2_s - w2_c));
            }
            else // om_s = 0 and w_c = 0
            {
                Psi2(i,j,k) = - dt;
            }

            // Psi3
            Complex Psi3;
            if (w_c != 0.)
            {
                Psi3 = - I * (theta3_c - theta_c) / (dt * w_c);
            }
            else // w_c = 0
            {
                Psi3 = 1._rt;
            }

            // Y1 (multiplies i*([k] \times J) in the update equation for <B>)
            if ((om_s != 0.) || (w_c != 0.))
            {
                Y1(i,j,k) = (1._rt - Psi1(i,j,k) - I * w_c * Psi2(i,j,k)) / (ep0 * (om2_s - w2_c));
            }
            else // om_s = 0 and w_c = 0
            {
                Y1(i,j,k) = 13._rt * dt2 / (24._rt * ep0);
            }

            // Y2 (multiplies rho_new in the update equation for <E>)
            if ((om_s != 0.) && (w_c != 0.))
            {
                Y2(i,j,k) = I * c2 * (ep0 * om2_s * Y1(i,j,k) - Psi3 + Psi1(i,j,k))
                            / (ep0 * om2_s * (theta2_c - 1._rt));
            }
            else if ((om_s != 0.) && (w_c == 0.))
            {
                Y2(i,j,k) = I * c2 * (C1 - C3 - dt2 * om2_s) / (ep0 * dt2 * om4_s);
            }
            else if ((om_s == 0.) && (w_c != 0.))
            {
                Y2(i,j,k) = c2 * (9._rt * dt2 * w2_c * theta3_c - dt2 * w2_c * theta_c
                            - 24._rt * theta3_c + 24._rt * theta_c + I * 8._rt * dt * w_c
                            + I * 24._rt * dt * w_c * theta3_c - I * 8._rt * dt * w_c * theta_c)
                            / (8._rt * ep0 * dt * w3_c * (1._rt - theta2_c));
            }
            else // om_s = 0 and w_c = 0
            {
                Y2(i,j,k) = - I * 5._rt * c2 * dt2 / (24._rt * ep0);
            }

            // Y3 (multiplies rho_old in the update equation for <E>)
            if ((om_s != 0.) && (w_c != 0.))
            {
                Y3(i,j,k) = I * c2 * (Psi3 - Psi1(i,j,k) - ep0 * theta2_c * om2_s * Y1(i,j,k))
                            / (ep0 * om2_s * (theta2_c - 1._rt));
            }
            else if ((om_s != 0.) && (w_c == 0.))
            {
                Y3(i,j,k) = I * c2 * (C3 - C1 + dt * om2_s * (S3_om - S1_om)) / (ep0 * dt2 * om4_s);
            }
            else if ((om_s == 0.) && (w_c != 0.))
            {
                Y3(i,j,k) = c2 * (9._rt * dt2 * w2_c * theta3_c - dt2 * w2_c * theta_c
                            - 16._rt * theta5_c + 8._rt * theta3_c + 8._rt * theta_c
                            + I * 12._rt * dt * w_c * theta5_c + I * 8._rt * dt * w_c * theta3_c
                            - I * 4._rt * dt * w_c * theta_c + I * 8._rt * dt * w_c * theta2_c)
                            / (8._rt * ep0 * dt * w3_c * (theta2_c - 1._rt));
            }
            else // om_s = 0 and w_c = 0
            {
                Y3(i,j,k) = - I * c2 * dt2 / (3._rt * ep0);
            }

            // Y4 (multiplies J in the update equation for <E>)
            Y4(i,j,k) = (Psi2(i,j,k) + I * ep0 * w_c * Y1(i,j,k)) / ep0;
        });
    }
}

void
PsatdAlgorithm::CurrentCorrection (
    const int lev,
    SpectralFieldData& field_data,
    std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
    const std::unique_ptr<amrex::MultiFab>& rho)
{
    // Profiling
    BL_PROFILE("PsatdAlgorithm::CurrentCorrection");

    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of J and rho
    field_data.ForwardTransform(lev, *current[0], Idx::Jx, 0);
    field_data.ForwardTransform(lev, *current[1], Idx::Jy, 0);
    field_data.ForwardTransform(lev, *current[2], Idx::Jz, 0);
    field_data.ForwardTransform(lev, *rho, Idx::rho_old, 0);
    field_data.ForwardTransform(lev, *rho, Idx::rho_new, 1);

    // Loop over boxes
    for (amrex::MFIter mfi(field_data.fields); mfi.isValid(); ++mfi){

        const amrex::Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> fields = field_data.fields[mfi].array();

        // Extract pointers for the k vectors
        const amrex::Real* const modified_kx_arr = modified_kx_vec[mfi].dataPtr();
        const amrex::Real* const modified_kx_arr_c = modified_kx_vec_centered[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
        const amrex::Real* const modified_ky_arr = modified_ky_vec[mfi].dataPtr();
        const amrex::Real* const modified_ky_arr_c = modified_ky_vec_centered[mfi].dataPtr();
#endif
        const amrex::Real* const modified_kz_arr = modified_kz_vec[mfi].dataPtr();
        const amrex::Real* const modified_kz_arr_c = modified_kz_vec_centered[mfi].dataPtr();

        // Local copy of member variables before GPU loop
        const amrex::Real dt = m_dt;

        // Galilean velocity
        const amrex::Real vgx = m_v_galilean[0];
        const amrex::Real vgy = m_v_galilean[1];
        const amrex::Real vgz = m_v_galilean[2];

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Shortcuts for the values of J and rho
            const Complex Jx = fields(i,j,k,Idx::Jx);
            const Complex Jy = fields(i,j,k,Idx::Jy);
            const Complex Jz = fields(i,j,k,Idx::Jz);
            const Complex rho_old = fields(i,j,k,Idx::rho_old);
            const Complex rho_new = fields(i,j,k,Idx::rho_new);

            // k vector values, and coefficients
            const amrex::Real kx = modified_kx_arr[i];
            const amrex::Real kx_c = modified_kx_arr_c[i];
#if (AMREX_SPACEDIM == 3)
            const amrex::Real ky = modified_ky_arr[j];
            const amrex::Real kz = modified_kz_arr[k];
            const amrex::Real ky_c = modified_ky_arr_c[j];
            const amrex::Real kz_c = modified_kz_arr_c[k];
#else
            constexpr amrex::Real ky = 0._rt;
            const     amrex::Real kz = modified_kz_arr[j];
            constexpr amrex::Real ky_c = 0._rt;
            const     amrex::Real kz_c = modified_kz_arr_c[j];
#endif
            constexpr Complex I = Complex{0._rt, 1._rt};

            const amrex::Real k_norm = std::sqrt(kx * kx + ky * ky + kz * kz);

            // Correct J
            if (k_norm != 0._rt)
            {
                const Complex k_dot_J = kx * Jx + ky * Jy + kz * Jz;
                const amrex::Real k_dot_vg = kx_c * vgx + ky_c * vgy + kz_c * vgz;

                // k_dot_vg = 0 always with standard PSATD (zero Galilean velocity)
                if ( k_dot_vg != 0._rt )
                {
                    const Complex rho_old_mod = rho_old * amrex::exp(I * k_dot_vg * dt);
                    const Complex den = 1._rt - amrex::exp(I * k_dot_vg * dt);

                    fields(i,j,k,Idx::Jx) = Jx - (k_dot_J - k_dot_vg * (rho_new - rho_old_mod) / den)
                        * kx / (k_norm * k_norm);

                    fields(i,j,k,Idx::Jy) = Jy - (k_dot_J - k_dot_vg * (rho_new - rho_old_mod) / den)
                        * ky / (k_norm * k_norm);

                    fields(i,j,k,Idx::Jz) = Jz - (k_dot_J - k_dot_vg * (rho_new - rho_old_mod) / den)
                        * kz / (k_norm * k_norm);
                }

                else
                {
                    fields(i,j,k,Idx::Jx) = Jx - (k_dot_J - I * (rho_new - rho_old) / dt)
                        * kx / (k_norm * k_norm);

                    fields(i,j,k,Idx::Jy) = Jy - (k_dot_J - I * (rho_new - rho_old) / dt)
                        * ky / (k_norm * k_norm);

                    fields(i,j,k,Idx::Jz) = Jz - (k_dot_J - I * (rho_new - rho_old) / dt)
                        * kz / (k_norm * k_norm);
                }
            }
        });
    }

    // Backward Fourier transform of J
    field_data.BackwardTransform(lev, *current[0], Idx::Jx, 0);
    field_data.BackwardTransform(lev, *current[1], Idx::Jy, 0);
    field_data.BackwardTransform(lev, *current[2], Idx::Jz, 0);
}

void
PsatdAlgorithm::VayDeposition (
    const int lev,
    SpectralFieldData& field_data,
    std::array<std::unique_ptr<amrex::MultiFab>,3>& current)
{
    // Profiling
    BL_PROFILE("PsatdAlgorithm::VayDeposition()");

    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of D (temporarily stored in current):
    // D is nodal and does not match the staggering of J, therefore we pass the
    // actual staggering of D (IntVect(1)) to the ForwardTransform function
    field_data.ForwardTransform(lev, *current[0], Idx::Jx, 0, IntVect(1));
    field_data.ForwardTransform(lev, *current[1], Idx::Jy, 0, IntVect(1));
    field_data.ForwardTransform(lev, *current[2], Idx::Jz, 0, IntVect(1));

    // Loop over boxes
    for (amrex::MFIter mfi(field_data.fields); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> fields = field_data.fields[mfi].array();

        // Extract pointers for the modified k vectors
        const amrex::Real* const modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
        const amrex::Real* const modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* const modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Shortcuts for the values of D
            const Complex Dx = fields(i,j,k,Idx::Jx);
#if (AMREX_SPACEDIM == 3)
            const Complex Dy = fields(i,j,k,Idx::Jy);
#endif
            const Complex Dz = fields(i,j,k,Idx::Jz);

            // Imaginary unit
            constexpr Complex I = Complex{0._rt, 1._rt};

            // Modified k vector values
            const amrex::Real kx_mod = modified_kx_arr[i];
#if (AMREX_SPACEDIM == 3)
            const amrex::Real ky_mod = modified_ky_arr[j];
            const amrex::Real kz_mod = modified_kz_arr[k];
#else
            const amrex::Real kz_mod = modified_kz_arr[j];
#endif

            // Compute Jx
            if (kx_mod != 0._rt) fields(i,j,k,Idx::Jx) = I * Dx / kx_mod;
            else                 fields(i,j,k,Idx::Jx) = 0._rt;

#if (AMREX_SPACEDIM == 3)
            // Compute Jy
            if (ky_mod != 0._rt) fields(i,j,k,Idx::Jy) = I * Dy / ky_mod;
            else                 fields(i,j,k,Idx::Jy) = 0._rt;
#endif

            // Compute Jz
            if (kz_mod != 0._rt) fields(i,j,k,Idx::Jz) = I * Dz / kz_mod;
            else                 fields(i,j,k,Idx::Jz) = 0._rt;
        });
    }

    // Backward Fourier transform of J
    field_data.BackwardTransform(lev, *current[0], Idx::Jx, 0);
    field_data.BackwardTransform(lev, *current[1], Idx::Jy, 0);
    field_data.BackwardTransform(lev, *current[2], Idx::Jz, 0);
}

#endif // WARPX_USE_PSATD
