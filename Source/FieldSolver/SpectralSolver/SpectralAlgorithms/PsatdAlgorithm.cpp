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
    const Array<Real,3>& v_galilean,
    const Real dt,
    const bool update_with_rho,
    const bool time_averaging)
    // Initializer list
    : SpectralBaseAlgorithm(spectral_kspace, dm, norder_x, norder_y, norder_z, nodal),
    // Initialize the centered finite-order modified k vectors: these are computed
    // always with the assumption of centered grids (argument nodal = true),
    // for both nodal and staggered simulations
    modified_kx_vec_centered(spectral_kspace.getModifiedKComponent(dm, 0, norder_x, true)),
#if (AMREX_SPACEDIM==3)
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
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

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

    // Allocate these coefficients only with averaged Galilean PSATD
    if (time_averaging)
    {
        C1_coef = SpectralRealCoefficients(ba, dm, 1, 0);
        S1_coef = SpectralRealCoefficients(ba, dm, 1, 0);
        C3_coef = SpectralRealCoefficients(ba, dm, 1, 0);
        S3_coef = SpectralRealCoefficients(ba, dm, 1, 0);
        Psi1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Psi2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Psi3_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        A1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        A2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Rhoold_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Rhonew_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
        Jcoef_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    }

    // Initialize coefficients
    InitializeSpectralCoefficients(spectral_kspace, dm, dt);
}

void
PsatdAlgorithm::pushSpectralFields (SpectralFieldData& f) const
{
    const bool update_with_rho = m_update_with_rho;
    const bool time_averaging  = m_time_averaging;
    const bool is_galilean     = m_is_galilean;

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi)
    {
        const Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();

        // These coefficients are always allocated
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
        Array4<const Complex> X1_arr = X1_coef[mfi].array();
        Array4<const Complex> X2_arr = X2_coef[mfi].array();
        Array4<const Complex> X3_arr = X3_coef[mfi].array();

        Array4<const Complex> X4_arr;
        Array4<const Complex> T2_arr;
        if (is_galilean)
        {
            X4_arr = X4_coef[mfi].array();
            T2_arr = T2_coef[mfi].array();
        }

        // These coefficients are allocated only with averaged Galilean PSATD
        Array4<const Complex> Psi1_arr;
        Array4<const Complex> Psi2_arr;
        Array4<const Complex> A1_arr;
        Array4<const Complex> Rhonew_arr;
        Array4<const Complex> Rhoold_arr;
        Array4<const Complex> Jcoef_arr;

        if (time_averaging)
        {
            Psi1_arr = Psi1_coef[mfi].array();
            Psi2_arr = Psi2_coef[mfi].array();
            A1_arr = A1_coef[mfi].array();
            Rhonew_arr = Rhonew_coef[mfi].array();
            Rhoold_arr = Rhoold_coef[mfi].array();
            Jcoef_arr = Jcoef_coef[mfi].array();
        }

        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

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
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0._rt;
            const     Real kz = modified_kz_arr[j];
#endif

            // Physical constants and imaginary unit
            constexpr Real c2 = PhysConst::c * PhysConst::c;
            constexpr Real inv_ep0 = 1._rt / PhysConst::ep0;
            constexpr Complex I = Complex{0._rt, 1._rt};

            // These coefficients are initialized in the function InitializeSpectralCoefficients
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Complex X1 = X1_arr(i,j,k);
            const Complex X2 = X2_arr(i,j,k);
            const Complex X3 = X3_arr(i,j,k);
            const Complex X4 = (is_galilean) ? X4_arr(i,j,k) : - S_ck * inv_ep0;
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

            fields(i,j,k,Idx::Bx) = T2 * C * Bx_old - I * T2 * S_ck * (ky * Ez_old - kz * Ey_old)
                + I * X1 * (ky * Jz - kz * Jy);

            fields(i,j,k,Idx::By) = T2 * C * By_old - I * T2 * S_ck * (kz * Ex_old - kx * Ez_old)
                + I * X1 * (kz * Jx - kx * Jz);

            fields(i,j,k,Idx::Bz) = T2 * C * Bz_old - I * T2 * S_ck * (kx * Ey_old - ky * Ex_old)
                + I * X1 * (kx * Jy - ky * Jx);

            // Additional update equations for averaged Galilean algorithm

            if (time_averaging)
            {
                // These coefficients are initialized in the function InitializeSpectralCoefficients below
                const Complex Psi1 = Psi1_arr(i,j,k);
                const Complex Psi2 = Psi2_arr(i,j,k);
                const Complex A1 = A1_arr(i,j,k);
                const Complex CRhoold = Rhoold_arr(i,j,k);
                const Complex CRhonew = Rhonew_arr(i,j,k);
                const Complex Jcoef = Jcoef_arr(i,j,k);

                fields(i,j,k,AvgIdx::Ex_avg) = Psi1 * Ex_old
                    - I * c2 * Psi2 * (ky * Bz_old - kz * By_old)
                    + Jcoef * Jx + (CRhonew * rho_new + CRhoold * rho_old) * kx;

                fields(i,j,k,AvgIdx::Ey_avg) = Psi1 * Ey_old
                    - I * c2 * Psi2 * (kz * Bx_old - kx * Bz_old)
                    + Jcoef * Jy + (CRhonew * rho_new + CRhoold * rho_old) * ky;

                fields(i,j,k,AvgIdx::Ez_avg) = Psi1 * Ez_old
                    - I * c2 * Psi2 * (kx * By_old - ky * Bx_old)
                    + Jcoef * Jz + (CRhonew * rho_new + CRhoold * rho_old) * kz;

                fields(i,j,k,AvgIdx::Bx_avg) = Psi1 * Bx_old + I * Psi2 * (ky * Ez_old - kz * Ey_old)
                    + I * inv_ep0 * A1 * (ky * Jz - kz * Jy);

                fields(i,j,k,AvgIdx::By_avg) = Psi1 * By_old + I * Psi2 * (kz * Ex_old - kx * Ez_old)
                    + I * inv_ep0 * A1 * (kz * Jx - kx * Jz);

                fields(i,j,k,AvgIdx::Bz_avg) = Psi1 * Bz_old + I * Psi2 * (kx * Ey_old - ky * Ex_old)
                    + I * inv_ep0 * A1 * (kx * Jy - ky * Jx);
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
    const bool time_averaging  = m_time_averaging;
    const bool is_galilean     = m_is_galilean;

    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Loop over boxes and allocate the corresponding coefficients for each box
    for (MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
        const Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const Real* kx = modified_kx_vec[mfi].dataPtr();
        const Real* kx_c = modified_kx_vec_centered[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* ky = modified_ky_vec[mfi].dataPtr();
        const Real* ky_c = modified_ky_vec_centered[mfi].dataPtr();
#endif
        const Real* kz = modified_kz_vec[mfi].dataPtr();
        const Real* kz_c = modified_kz_vec_centered[mfi].dataPtr();

        // Coefficients always allocated
        Array4<Real> C = C_coef[mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();
        Array4<Complex> X1 = X1_coef[mfi].array();
        Array4<Complex> X2 = X2_coef[mfi].array();
        Array4<Complex> X3 = X3_coef[mfi].array();

        Array4<Complex> X4;
        Array4<Complex> T2;
        if (is_galilean)
        {
            X4 = X4_coef[mfi].array();
            T2 = T2_coef[mfi].array();
        }

        // Coefficients allocated only with averaged Galilean PSATD
        Array4<Real> C1;
        Array4<Real> S1;
        Array4<Real> C3;
        Array4<Real> S3;
        Array4<Complex> Psi1;
        Array4<Complex> Psi2;
        Array4<Complex> Psi3;
        Array4<Complex> A1;
        Array4<Complex> A2;
        Array4<Complex> CRhoold;
        Array4<Complex> CRhonew;
        Array4<Complex> Jcoef;

        if (time_averaging)
        {
            C1 = C1_coef[mfi].array();
            S1 = S1_coef[mfi].array();
            C3 = C3_coef[mfi].array();
            S3 = S3_coef[mfi].array();
            Psi1 = Psi1_coef[mfi].array();
            Psi2 = Psi2_coef[mfi].array();
            Psi3 = Psi3_coef[mfi].array();
            A1 = A1_coef[mfi].array();
            A2 = A2_coef[mfi].array();
            CRhoold = Rhoold_coef[mfi].array();
            CRhonew = Rhonew_coef[mfi].array();
            Jcoef = Jcoef_coef[mfi].array();
        }

        // Extract Galilean velocity
        Real vx = m_v_galilean[0];
#if (AMREX_SPACEDIM==3)
        Real vy = m_v_galilean[1];
#endif
        Real vz = m_v_galilean[2];

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of k vector
            const Real knorm = std::sqrt(
                std::pow(kx[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(ky[j], 2) + std::pow(kz[k], 2));
#else
                std::pow(kz[j], 2));
#endif
            // Calculate norm of k vector (centered)
            const Real knorm_c = std::sqrt(
                std::pow(kx_c[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(ky_c[j], 2) + std::pow(kz_c[k], 2));
#else
                std::pow(kz_c[j], 2));
#endif
            // Physical constants and imaginary unit
            constexpr Real c = PhysConst::c;
            constexpr Real c2 = c*c;
            constexpr Real ep0 = PhysConst::ep0;
            constexpr Complex I = Complex{0._rt, 1._rt};

            // Auxiliary coefficients used when update_with_rho=false
            const Real dt2 = dt * dt;
            const Real dt3 = dt * dt2;
            Complex X2_old, X3_old;

            // Calculate the dot product of the k vector with the Galilean velocity.
            // This has to be computed always with the centered (that is, nodal) finite-order
            // modified k vectors, to work correctly for both nodal and staggered simulations.
            // kv = 0 always with standard PSATD (zero Galilean velocity).
            const Real kv = kx_c[i]*vx +
#if (AMREX_SPACEDIM==3)
                ky_c[j]*vy + kz_c[k]*vz;
#else
                kz_c[j]*vz;
#endif

            // Note that:
            // - X1 multiplies i*(k \times J) in the update equation for B
            // - X2 multiplies rho_new if update_with_rho = 1 or (k \dot E)
            //      if update_with_rho = 0 in the update equation for E
            // - X3 multiplies rho_old if update_with_rho = 1 or (k \dot J)
            //      if update_with_rho = 0 in the update equation for E
            // - X4 multiplies J in the update equation for E

            if (knorm != 0. && knorm_c != 0.)
            {
                // Auxiliary coefficients
                const Real om = c * knorm;
                const Real om2 = om * om;
                const Real om3 = om * om2;
                const Real om_c = c * knorm_c;
                const Real om2_c = om_c * om_c;
                const Real om3_c = om_c * om2_c;
                const Complex tmp1 = amrex::exp(  I * om * dt);
                const Complex tmp2 = amrex::exp(- I * om * dt);

                C   (i,j,k) = std::cos(om * dt);
                S_ck(i,j,k) = std::sin(om * dt) / om;

                if (time_averaging)
                {
                    C1(i,j,k) = std::cos(0.5_rt * om * dt);
                    S1(i,j,k) = std::sin(0.5_rt * om * dt);
                    C3(i,j,k) = std::cos(1.5_rt * om * dt);
                    S3(i,j,k) = std::sin(1.5_rt * om * dt);
                }

                const Real nu = kv / om_c;
                const Real nu2 = nu * nu;
                const Complex theta = amrex::exp(I * nu * om_c * dt * 0.5_rt);
                const Complex theta_star = amrex::exp(- I * nu * om_c * dt * 0.5_rt);

                // T2 = 1 always with standard PSATD (zero Galilean velocity)
                if (is_galilean)
                {
                    T2(i,j,k) = theta * theta;
                }
                const Complex T2_tmp = (is_galilean) ? T2(i,j,k) : 1.0_rt;

                // nu = 0 always with standard PSATD (zero Galilean velocity): skip this block
                if (nu != om/om_c && nu != -om/om_c && nu != 0.)
                {
                    // x1 is the coefficient chi_1 in equation (12c)
                    Complex x1 = om2_c / (om2 - nu2 * om2_c)
                        * (theta_star - theta * C(i,j,k) + I * nu * om_c * theta * S_ck(i,j,k));

                    X1(i,j,k) = theta * x1 / (ep0 * om2_c);

                    if (update_with_rho)
                    {
                        X2(i,j,k) = c2 * (x1 * om2 - theta * (1._rt - C(i,j,k)) * om2_c)
                            / (theta_star - theta) / (ep0 * om2_c * om2);

                        X3(i,j,k) = c2 * (x1 * om2 - theta_star * (1._rt - C(i,j,k)) * om2_c)
                            / (theta_star - theta) / (ep0 * om2_c * om2);
                    }

                    else // update_with_rho = 0
                    {
                        X2_old = (x1 * om2 - theta * (1._rt - C(i,j,k)) * om2_c)
                            / (theta_star - theta);

                        X3_old = (x1 * om2 - theta_star * (1._rt - C(i,j,k)) * om2_c)
                            / (theta_star - theta);

                        X2(i,j,k) = c2 * T2_tmp * (X2_old - X3_old) / (om2_c * om2);

                        X3(i,j,k) = I * c2 * X2_old * (T2_tmp - 1._rt) / (ep0 * nu * om3_c * om2);
                    }

                    if (is_galilean)
                    {
                        X4(i,j,k) = I * nu * om_c * X1(i,j,k) - T2_tmp * S_ck(i,j,k) / ep0;
                    }

                    // Averaged Galilean algorithm
                    if (time_averaging)
                    {
                        Complex C_rho = I * c2 / ((1._rt - T2_tmp) * ep0);

                        Psi1(i,j,k) = theta * ((om * S1(i,j,k) + I * nu * om_c * C1(i,j,k))
                            - T2_tmp * (om * S3(i,j,k) + I * nu * om_c * C3(i,j,k)))
                            / (dt * (nu2 * om2_c - om2));

                        Psi2(i,j,k) = theta * ((om * C1(i,j,k) - I * nu * om_c * S1(i,j,k))
                            - T2_tmp * (om * C3(i,j,k) - I * nu * om_c * S3(i,j,k)))
                            / (om * dt * (nu2 * om2_c - om2));

                        Psi3(i,j,k) = I * theta * (1._rt - T2_tmp) / (nu * om_c * dt);

                        A1(i,j,k) = (Psi1(i,j,k) - 1._rt + I * nu * om_c * Psi2(i,j,k))
                            / (nu2 * om2_c - om2);

                        A2(i,j,k) = (Psi3(i,j,k) - Psi1(i,j,k)) / om2;

                        CRhoold(i,j,k) = C_rho * (T2_tmp * A1(i,j,k) - A2(i,j,k));

                        CRhonew(i,j,k) = C_rho * (A2(i,j,k) - A1(i,j,k));

                        Jcoef(i,j,k) = (I * nu * om_c * A1(i,j,k) + Psi2(i,j,k)) / ep0;
                    }
                }

                // nu = 0 always with standard PSATD (zero Galilean velocity)
                if (nu == 0.)
                {
                    X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0 * om2);

                    if (update_with_rho)
                    {
                        X2(i,j,k) = c2 * (1._rt - S_ck(i,j,k) / dt) / (ep0 * om2);

                        X3(i,j,k) = c2 * (C(i,j,k) - S_ck(i,j,k) / dt) / (ep0 * om2);
                    }

                    else // update_with_rho = 0
                    {
                        X2(i,j,k) = c2 * (1._rt - C(i,j,k)) / om2;

                        X3(i,j,k) = c2 * (S_ck(i,j,k) / dt - 1._rt) * dt / (ep0 * om2);
                    }

                    if (is_galilean)
                    {
                        X4(i,j,k) = - S_ck(i,j,k) / ep0;
                    }

                    // Averaged Galilean algorithm
                    if (time_averaging)
                    {
                        Psi1(i,j,k) = (S3(i,j,k) - S1(i,j,k)) / (om * dt);

                        Psi2(i,j,k) = (C3(i,j,k) - C1(i,j,k)) / (om2 * dt);

                        Psi3(i,j,k) = 1._rt;

                        A1(i,j,k) = (om * dt + S1(i,j,k) - S3(i,j,k)) / (om3 * dt);

                        A2(i,j,k) = (om * dt + S1(i,j,k) - S3(i,j,k)) / (om3 * dt);

                        CRhoold(i,j,k) = 2._rt * I * c2 * S1(i,j,k) * (dt * C(i,j,k) - S_ck(i,j,k))
                            / (om3 * dt2 * ep0);

                        CRhonew(i,j,k) = - I * c2 * (om2 * dt2 - C1(i,j,k) + C3(i,j,k))
                            / (om2 * om2 * dt2 * ep0);

                        Jcoef(i,j,k) = (I * nu * om_c * A1(i,j,k) + Psi2(i,j,k)) / ep0;
                    }
                }

                // nu = 0 always with standard PSATD (zero Galilean velocity): skip this block
                if (nu == om/om_c)
                {
                    X1(i,j,k) = (1._rt - tmp1 * tmp1 + 2._rt * I * om * dt) / (4._rt * ep0 * om2);

                    if (update_with_rho)
                    {
                        X2(i,j,k) = c2 * (- 3._rt + 4._rt * tmp1 - tmp1 * tmp1 - 2._rt * I * om * dt)
                            / (4._rt * ep0 * om2 * (tmp1 - 1._rt));

                        X3(i,j,k) = c2 * (3._rt - 2._rt * tmp2 - 2._rt * tmp1 + tmp1 * tmp1
                            - 2._rt * I * om * dt) / (4._rt * ep0 * om2 * (tmp1 - 1._rt));
                    }

                    else // update_with_rho = 0
                    {
                        X2(i,j,k) = c2 * (1._rt - C(i,j,k)) * tmp1 / om2;

                        X3(i,j,k) = c2 * (2._rt * om * dt - I * tmp1 * tmp1 + 4._rt * I * tmp1 - 3._rt * I)
                            / (4._rt * ep0 * om3);
                    }

                    if (is_galilean)
                    {
                        X4(i,j,k) = (- I + I * tmp1 * tmp1 - 2._rt * om * dt) / (4._rt * ep0 * om);
                    }

                    // Averaged Galilean algorithm
                    if (time_averaging)
                    {
                        Complex C_rho = I * c2 / ((1._rt - T2_tmp) * ep0);

                        Psi1(i,j,k) = (2._rt * om * dt + I * tmp1 - I * tmp1 * tmp1 * tmp1)
                            / (4._rt * om * dt);

                        Psi2(i,j,k) = (- 2._rt * I * om * dt - tmp1 + tmp1 * tmp1 * tmp1)
                            / (4._rt * om2 * dt);

                        Psi3(i,j,k) = I * theta * (1._rt - T2_tmp) / (nu * om_c * dt);

                        A1(i,j,k) = (2._rt * om * dt + I * (4._rt * om2 * dt2 - tmp1 + tmp1 * tmp1 * tmp1))
                            / (8._rt * om3 * dt);

                        A2(i,j,k) = (Psi3(i,j,k) - Psi1(i,j,k)) / om2;

                        CRhoold(i,j,k) = C_rho * (T2_tmp * A1(i,j,k) - A2(i,j,k));

                        CRhonew(i,j,k) = C_rho * (A2(i,j,k) - A1(i,j,k));

                        Jcoef(i,j,k) = (I * nu * om_c * A1(i,j,k) + Psi2(i,j,k)) / ep0;
                    }
                }

                // nu = 0 always with standard PSATD (zero Galilean velocity): skip this block
                if (nu == -om/om_c)
                {
                    X1(i,j,k) = (1._rt - tmp2 * tmp2 - 2._rt * I * om * dt) / (4._rt * ep0 * om2);

                    if (update_with_rho)
                    {
                        X2(i,j,k) = c2 * (- 4._rt + 3._rt * tmp1 + tmp2 - 2._rt * I * om * dt * tmp1)
                            / (4._rt * ep0 * om2 * (tmp1 - 1._rt));

                        X3(i,j,k) = c2 * (2._rt - tmp2 - 3._rt * tmp1 + 2._rt * tmp1 * tmp1
                            - 2._rt * I * om * dt * tmp1) / (4._rt * ep0 * om2 * (tmp1 - 1._rt));
                    }

                    else // update_with_rho = 0
                    {
                        X2(i,j,k) = c2 * (1._rt - C(i,j,k)) * tmp2 / om2;

                        X3(i,j,k) = c2 * (2._rt * om * dt + I * tmp2 * tmp2 - 4._rt * I * tmp2 + 3._rt * I)
                            / (4._rt * ep0 * om3);
                    }

                    if (is_galilean)
                    {
                        X4(i,j,k) = (I - I * tmp2 * tmp2 - 2._rt * om * dt) / (4._rt * ep0 * om);
                    }

                    // Averaged Galilean algorithm
                    if (time_averaging)
                    {
                        Complex C_rho = I * c2 / ((1._rt - T2_tmp) * ep0);

                        Psi1(i,j,k) = (2._rt * om * dt - I * tmp2 + I * tmp2 * tmp2 * tmp2)
                            / (4._rt * om * dt);

                        Psi2(i,j,k) = (2._rt * I * om * dt - tmp2 + tmp2 * tmp2 * tmp2)
                            / (4._rt * om2 * dt);

                        Psi3(i,j,k) = I * theta * (1._rt - T2_tmp) / (nu * om_c * dt);

                        A1(i,j,k) = (2._rt * om * dt * (1._rt - 2._rt * I * om * dt)
                            + I * (tmp2 - tmp2 * tmp2 * tmp2)) / (8._rt * om3 * dt);

                        A2(i,j,k) = (Psi3(i,j,k) - Psi1(i,j,k)) / om2;

                        CRhoold(i,j,k) = C_rho * (T2_tmp * A1(i,j,k) - A2(i,j,k));

                        CRhonew(i,j,k) = C_rho * (A2(i,j,k) - A1(i,j,k));

                        Jcoef(i,j,k) = (I * nu * om_c * A1(i,j,k) + Psi2(i,j,k)) / ep0;
                    }
                }
            }

            else if (knorm != 0. && knorm_c == 0.)
            {
                const Real om = c * knorm;
                const Real om2 = om * om;
                const Real om3 = om * om2;

                C(i,j,k) = std::cos(om * dt);

                S_ck(i,j,k) = std::sin(om * dt) / om;

                X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0 * om2);

                if (update_with_rho)
                {
                    X2(i,j,k) = c2 * (1._rt - S_ck(i,j,k) / dt) / (ep0 * om2);

                    X3(i,j,k) = c2 * (C(i,j,k) - S_ck(i,j,k) / dt) / (ep0 * om2);
                }

                else // update_with_rho = 0
                {
                    X2(i,j,k) = c2 * (1._rt - C(i,j,k)) / om2;

                    X3(i,j,k) = c2 * (S_ck(i,j,k) / dt - 1._rt) * dt / (ep0 * om2);
                }

                if (is_galilean)
                {
                    X4(i,j,k) = - S_ck(i,j,k) / ep0;

                    T2(i,j,k) = 1._rt;
                }

                // Averaged Galilean algorithm
                if (time_averaging)
                {
                    C1(i,j,k) = std::cos(0.5_rt * om * dt);
                    S1(i,j,k) = std::sin(0.5_rt * om * dt);
                    C3(i,j,k) = std::cos(1.5_rt * om * dt);
                    S3(i,j,k) = std::sin(1.5_rt * om * dt);

                    Psi1(i,j,k) = (S3(i,j,k) - S1(i,j,k)) / (om * dt);

                    Psi2(i,j,k) = (C3(i,j,k) - C1(i,j,k)) / (om2 * dt);

                    Psi3(i,j,k) = 1._rt;

                    A1(i,j,k) = (om * dt + S1(i,j,k) - S3(i,j,k)) / (om3 * dt);

                    A2(i,j,k) = (om * dt + S1(i,j,k) - S3(i,j,k)) / (om3 * dt);

                    CRhoold(i,j,k) = 2._rt * I * c2 * S1(i,j,k) * (dt * C(i,j,k) - S_ck(i,j,k))
                        / (om3 * dt2 * ep0);

                    CRhonew(i,j,k) = - I * c2 * (om2 * dt2 - C1(i,j,k) + C3(i,j,k))
                        / (om2 * om2 * dt2 * ep0);

                    Jcoef(i,j,k) = Psi2(i,j,k) / ep0;
                }
            }

            else if (knorm == 0. && knorm_c != 0.)
            {
                const Real om_c = c * knorm_c;
                const Real om2_c = om_c * om_c;
                const Real om3_c = om_c * om2_c;
                const Real nu  = kv / om_c;
                const Real nu2 = nu * nu;
                const Complex theta = amrex::exp(I * nu * om_c * dt * 0.5_rt);

                C(i,j,k) = 1._rt;

                S_ck(i,j,k) = dt;

                if (is_galilean)
                {
                    T2(i,j,k) = theta * theta;
                }
                const Complex T2_tmp = (is_galilean) ? T2(i,j,k) : 1.0_rt;

                // nu = 0 always with standard PSATD (zero Galilean velocity): skip this block
                if (nu != 0.)
                {
                    X1(i,j,k) = (- 1._rt + T2_tmp - I * nu * om_c * dt * T2_tmp) / (ep0 * nu2 * om2_c);

                    if (update_with_rho)
                    {
                        X2(i,j,k) = c2 * (1._rt - T2_tmp + I * nu * om_c * dt * T2_tmp
                            + 0.5_rt * nu2 * om2_c * dt2 * T2_tmp) / (ep0 * nu2 * om2_c * (T2_tmp - 1._rt));

                        X3(i,j,k) = c2 * (1._rt - T2_tmp + I * nu * om_c * dt * T2_tmp
                            + 0.5_rt * nu2 * om2_c * dt2) / (ep0 * nu2 * om2_c * (T2_tmp - 1._rt));
                    }

                    else // update_with_rho = 0
                    {
                        X2(i,j,k) = c2 * dt2 * T2_tmp * 0.5_rt;

                        X3(i,j,k) = c2 * (2._rt * I - 2._rt * nu * om_c * dt * T2_tmp
                            + I * nu2 * om2_c * dt2 * T2_tmp) / (2._rt * ep0 * nu2 * nu * om3_c);
                    }

                    if (is_galilean)
                    {
                        X4(i,j,k) = I * (T2_tmp - 1._rt) / (ep0 * nu * om_c);
                    }

                    // Averaged Galilean algorithm
                    if (time_averaging)
                    {
                        Complex C_rho = I * c2 / ((1._rt - T2_tmp) * ep0);

                        Psi1(i,j,k) = I * theta * (1._rt - T2_tmp) / (nu * om_c * dt);

                        Psi2(i,j,k) = theta * (2._rt - I * nu * om_c * dt + T2_tmp * (3._rt * I * nu * om_c * dt - 2._rt))
                            / (2._rt * nu2 * om2_c * dt);

                        Psi3(i,j,k) = I * theta * (1._rt - T2_tmp) / (nu * om_c * dt);

                        A1(i,j,k) = (Psi1(i,j,k) - 1._rt + I * nu * om_c * Psi2(i,j,k)) / (nu2 * om2_c);

                        A2(i,j,k) = theta * (8._rt * I * (T2_tmp - 1._rt) + 4._rt * nu * om_c * dt
                            * (3._rt * T2_tmp - 1._rt) + I * nu2 * om2_c * dt2 * (1._rt - 9._rt * T2_tmp))
                            / (8._rt * nu2 * nu * om2_c * om_c * dt);

                        CRhoold(i,j,k) = C_rho * (T2_tmp * A1(i,j,k) - A2(i,j,k));

                        CRhonew(i,j,k) = C_rho * (A2(i,j,k) - A1(i,j,k));

                        Jcoef(i,j,k) = (I * nu * om_c * A1(i,j,k) + Psi2(i,j,k)) / ep0;
                    }
                }

                else // nu = 0
                {
                    X1(i,j,k) = dt2 / (2._rt * ep0);

                    if (update_with_rho)
                    {
                        X2(i,j,k) = c2 * dt2 / (6._rt * ep0);

                        X3(i,j,k) = - c2 * dt2 / (3._rt * ep0);
                    }

                    else // update_with_rho = 0
                    {
                        X2(i,j,k) = c2 * dt2 * 0.5_rt;

                        X3(i,j,k) = - c2 * dt3 / (6._rt * ep0);
                    }

                    if (is_galilean)
                    {
                        X4(i,j,k) = - dt / ep0;
                    }

                    // Averaged Galilean algorithm
                    if (time_averaging)
                    {
                        Psi1(i,j,k) = 1._rt;

                        Psi2(i,j,k) = - dt;

                        Psi3(i,j,k) = 1._rt;

                        A1(i,j,k) = 13._rt * dt2 / 24._rt;

                        A2(i,j,k) = 13._rt * dt2 / 24._rt;

                        CRhoold(i,j,k) = - I * c2 * dt2 / (3._rt * ep0);

                        CRhonew(i,j,k) = - 5._rt * I * c2 * dt2 / (24._rt * ep0);

                        Jcoef(i,j,k) = - dt / ep0;
                    }
                }
            }

            else if (knorm == 0. && knorm_c == 0.)
            {
                C(i,j,k) = 1._rt;

                S_ck(i,j,k) = dt;

                X1(i,j,k) = dt2 / (2._rt * ep0);

                if (update_with_rho)
                {
                    X2(i,j,k) = c2 * dt2 / (6._rt * ep0);

                    X3(i,j,k) = - c2 * dt2 / (3._rt * ep0);
                }

                else // update_with_rho = 0
                {
                    X2(i,j,k) = c2 * dt2 * 0.5_rt;

                    X3(i,j,k) = - c2 * dt3 / (6._rt * ep0);
                }

                // T2 = 1 always with standard PSATD (zero Galilean velocity)
                if (is_galilean)
                {
                    X4(i,j,k) = - dt / ep0;

                    T2(i,j,k) = 1._rt;
                }

                // Averaged Galilean algorithm
                if (time_averaging)
                {
                    Psi1(i,j,k) = 1._rt;

                    Psi2(i,j,k) = - dt;

                    Psi3(i,j,k) = 1._rt;

                    A1(i,j,k) = 13._rt * dt2 / 24._rt;

                    A2(i,j,k) = 13._rt * dt2 / 24._rt;

                    CRhoold(i,j,k) = - I * c2 * dt2 / (3._rt * ep0);

                    CRhonew(i,j,k) = - 5._rt * I * c2 * dt2 / (24._rt * ep0);

                    Jcoef(i,j,k) = - dt / ep0;
                }
            }
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
#if (AMREX_SPACEDIM==3)
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
#if (AMREX_SPACEDIM==3)
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
#if (AMREX_SPACEDIM==3)
        const amrex::Real* const modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* const modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Shortcuts for the values of D
            const Complex Dx = fields(i,j,k,Idx::Jx);
#if (AMREX_SPACEDIM==3)
            const Complex Dy = fields(i,j,k,Idx::Jy);
#endif
            const Complex Dz = fields(i,j,k,Idx::Jz);

            // Imaginary unit
            constexpr Complex I = Complex{0._rt, 1._rt};

            // Modified k vector values
            const amrex::Real kx_mod = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const amrex::Real ky_mod = modified_ky_arr[j];
            const amrex::Real kz_mod = modified_kz_arr[k];
#else
            const amrex::Real kz_mod = modified_kz_arr[j];
#endif

            // Compute Jx
            if (kx_mod != 0._rt) fields(i,j,k,Idx::Jx) = I * Dx / kx_mod;
            else                 fields(i,j,k,Idx::Jx) = 0._rt;

#if (AMREX_SPACEDIM==3)
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
