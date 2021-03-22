#include "ComovingPsatdAlgorithm.H"
#include "Utils/WarpXConst.H"

#if WARPX_USE_PSATD

using namespace amrex;

ComovingPsatdAlgorithm::ComovingPsatdAlgorithm (const SpectralKSpace& spectral_kspace,
                                                const DistributionMapping& dm,
                                                const int norder_x, const int norder_y,
                                                const int norder_z, const bool nodal,
                                                const amrex::Array<amrex::Real, 3>& v_comoving,
                                                const amrex::Real dt,
                                                const bool update_with_rho)
     // Members initialization
     : SpectralBaseAlgorithm(spectral_kspace, dm, norder_x, norder_y, norder_z, nodal),
       // Initialize the infinite-order k vectors (the argument n_order = -1 selects
       // the infinite order option, the argument nodal = false is then irrelevant)
       kx_vec(spectral_kspace.getModifiedKComponent(dm, 0, -1, false)),
#if (AMREX_SPACEDIM==3)
       ky_vec(spectral_kspace.getModifiedKComponent(dm, 1, -1, false)),
       kz_vec(spectral_kspace.getModifiedKComponent(dm, 2, -1, false)),
#else
       kz_vec(spectral_kspace.getModifiedKComponent(dm, 1, -1, false)),
#endif
       m_v_comoving(v_comoving),
       m_dt(dt),
       m_update_with_rho(update_with_rho)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate arrays of real spectral coefficients
    C_coef    = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);

    // Allocate arrays of complex spectral coefficients
    X1_coef     = SpectralComplexCoefficients(ba, dm, 1, 0);
    X2_coef     = SpectralComplexCoefficients(ba, dm, 1, 0);
    X3_coef     = SpectralComplexCoefficients(ba, dm, 1, 0);
    X4_coef     = SpectralComplexCoefficients(ba, dm, 1, 0);
    Theta2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

    // Initialize real and complex spectral coefficients
    InitializeSpectralCoefficients(spectral_kspace, dm, dt);
}

void
ComovingPsatdAlgorithm::pushSpectralFields (SpectralFieldData& f) const
{
    // Loop over boxes
    for (amrex::MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        const amrex::Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> fields = f.fields[mfi].array();

        // Extract arrays for the coefficients
        amrex::Array4<const amrex::Real> C_arr    = C_coef   [mfi].array();
        amrex::Array4<const amrex::Real> S_ck_arr = S_ck_coef[mfi].array();
        amrex::Array4<const Complex>     X1_arr   = X1_coef  [mfi].array();
        amrex::Array4<const Complex>     X2_arr   = X2_coef  [mfi].array();
        amrex::Array4<const Complex>     X3_arr   = X3_coef  [mfi].array();
        amrex::Array4<const Complex>     X4_arr   = X4_coef  [mfi].array();

        // Extract pointers for the k vectors
        const amrex::Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const amrex::Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Record old values of the fields to be updated
            using Idx = SpectralFieldIndex;
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
            const amrex::Real kx_mod = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const amrex::Real ky_mod = modified_ky_arr[j];
            const amrex::Real kz_mod = modified_kz_arr[k];
#else
            constexpr amrex::Real ky_mod = 0._rt;
            const     amrex::Real kz_mod = modified_kz_arr[j];
#endif

            // Physical constant c**2 and imaginary unit
            constexpr amrex::Real c2 = PhysConst::c*PhysConst::c;
            constexpr Complex I = Complex{0._rt,1._rt};

            // The definition of these coefficients is explained in more detail
            // in the function InitializeSpectralCoefficients below
            const amrex::Real C    = C_arr(i,j,k);
            const amrex::Real S_ck = S_ck_arr(i,j,k);
            const Complex     X1   = X1_arr(i,j,k);
            const Complex     X2   = X2_arr(i,j,k);
            const Complex     X3   = X3_arr(i,j,k);
            const Complex     X4   = X4_arr(i,j,k);

            // Update E
            fields(i,j,k,Idx::Ex) = C*Ex_old + S_ck*c2*I*(ky_mod*Bz_old - kz_mod*By_old)
                + X4*Jx - I*(X2*rho_new - X3*rho_old)*kx_mod;

            fields(i,j,k,Idx::Ey) = C*Ey_old + S_ck*c2*I*(kz_mod*Bx_old - kx_mod*Bz_old)
                + X4*Jy - I*(X2*rho_new - X3*rho_old)*ky_mod;

            fields(i,j,k,Idx::Ez) = C*Ez_old + S_ck*c2*I*(kx_mod*By_old - ky_mod*Bx_old)
                + X4*Jz - I*(X2*rho_new - X3*rho_old)*kz_mod;

            // Update B
            fields(i,j,k,Idx::Bx) = C*Bx_old - S_ck*I*(ky_mod*Ez_old - kz_mod*Ey_old)
                + X1*I*(ky_mod*Jz - kz_mod*Jy);

            fields(i,j,k,Idx::By) = C*By_old - S_ck*I*(kz_mod*Ex_old - kx_mod*Ez_old)
                + X1*I*(kz_mod*Jx - kx_mod*Jz);

            fields(i,j,k,Idx::Bz) = C*Bz_old - S_ck*I*(kx_mod*Ey_old - ky_mod*Ex_old)
                + X1*I*(kx_mod*Jy - ky_mod*Jx);
        });
    }
}

void ComovingPsatdAlgorithm::InitializeSpectralCoefficients (const SpectralKSpace& spectral_kspace,
                                                             const amrex::DistributionMapping& dm,
                                                             const amrex::Real dt)
{
    const amrex::BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Loop over boxes and allocate the corresponding coefficients for each box
    for (amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const amrex::Real* kx_mod = modified_kx_vec[mfi].dataPtr();
        const amrex::Real* kx     = kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const amrex::Real* ky_mod = modified_ky_vec[mfi].dataPtr();
        const amrex::Real* ky     = ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* kz_mod = modified_kz_vec[mfi].dataPtr();
        const amrex::Real* kz     = kz_vec[mfi].dataPtr();

        // Extract arrays for the coefficients
        amrex::Array4<amrex::Real> C    = C_coef     [mfi].array();
        amrex::Array4<amrex::Real> S_ck = S_ck_coef  [mfi].array();
        amrex::Array4<Complex>     X1   = X1_coef    [mfi].array();
        amrex::Array4<Complex>     X2   = X2_coef    [mfi].array();
        amrex::Array4<Complex>     X3   = X3_coef    [mfi].array();
        amrex::Array4<Complex>     X4   = X4_coef    [mfi].array();
        amrex::Array4<Complex>     T2   = Theta2_coef[mfi].array();

        // Store comoving velocity
        const amrex::Real vx = m_v_comoving[0];
#if (AMREX_SPACEDIM==3)
        const amrex::Real vy = m_v_comoving[1];
#endif
        const amrex::Real vz = m_v_comoving[2];

        // Loop over indices within one box
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of finite-order k vector
            const amrex::Real knorm_mod = std::sqrt(
                std::pow(kx_mod[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(ky_mod[j], 2) +
                std::pow(kz_mod[k], 2));
#else
                std::pow(kz_mod[j], 2));
#endif
            // Calculate norm of infinite-order k vector
            const amrex::Real knorm = std::sqrt(
                std::pow(kx[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(ky[j], 2) +
                std::pow(kz[k], 2));
#else
                std::pow(kz[j], 2));
#endif
            // Physical constants c, c**2, and epsilon_0, and imaginary unit
            constexpr amrex::Real c   = PhysConst::c;
            constexpr amrex::Real c2  = c*c;
            constexpr amrex::Real ep0 = PhysConst::ep0;
            constexpr Complex     I   = Complex{0._rt, 1._rt};

            // Auxiliary coefficients used when update_with_rho=false
            const amrex::Real dt2 = dt * dt;

            // Calculate dot product of k vector with comoving velocity
            const amrex::Real kv = kx[i]*vx +
#if (AMREX_SPACEDIM==3)
                ky[j]*vy + kz[k]*vz;
#else
                kz[j]*vz;
#endif

            if (knorm_mod != 0. && knorm != 0.) {

                // Auxiliary coefficients
                const amrex::Real om_mod  = c * knorm_mod;
                const amrex::Real om2_mod = om_mod * om_mod;
                const amrex::Real om  = c * knorm;
                const amrex::Real om2 = om * om;
                const Complex tmp1 = amrex::exp(  I * om_mod * dt);
                const Complex tmp2 = amrex::exp(- I * om_mod * dt);
                const Complex tmp1_sqrt = amrex::exp(  I * om_mod * dt * 0.5_rt);
                const Complex tmp2_sqrt = amrex::exp(- I * om_mod * dt * 0.5_rt);

                C   (i,j,k) = std::cos(om_mod * dt);
                S_ck(i,j,k) = std::sin(om_mod * dt) / om_mod;

                const amrex::Real nu = - kv / om;
                const Complex theta      = amrex::exp(  I * nu * om * dt * 0.5_rt);
                const Complex theta_star = amrex::exp(- I * nu * om * dt * 0.5_rt);

                T2(i,j,k) = theta * theta;

                if ( (nu != om_mod/om) && (nu != -om_mod/om) && (nu != 0.) ) {

                    Complex x1 = om2 / (om2_mod - nu * nu * om2)
                        * (theta_star - theta * C(i,j,k) + I * nu * om * theta * S_ck(i,j,k));

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = x1 / (ep0 * om2);

                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * (x1 * om2_mod - theta * (1._rt - C(i,j,k)) * om2)
                        / (theta_star - theta) / (ep0 * om2 * om2_mod);
                    X3(i,j,k) = c2 * (x1 * om2_mod - theta_star * (1._rt - C(i,j,k)) * om2)
                        / (theta_star - theta) / (ep0 * om2 * om2_mod);

                    // X4 multiplies J in the update equation for E
                    X4(i,j,k) = I * nu * om * X1(i,j,k) - theta * S_ck(i,j,k) / ep0;
                }

                // Limits for nu = 0
                if (nu == 0.) {

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0 * om2_mod);

                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * (1._rt - S_ck(i,j,k) / dt) / (ep0 * om2_mod);
                    X3(i,j,k) = c2 * (C(i,j,k) - S_ck(i,j,k) / dt) / (ep0 * om2_mod);

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = - S_ck(i,j,k) / ep0;
                }

                // Limits for nu = omega_mod/omega
                if (nu == om_mod/om) {

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = tmp1_sqrt * (1._rt - tmp2 * tmp2 - 2._rt * I * om_mod * dt) / (4._rt * ep0 * om2_mod);

                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * (- 4._rt + 3._rt * tmp1 + tmp2 - 2._rt * I * om_mod * dt * tmp1)
                        / (4._rt * ep0 * om2_mod * (tmp1 - 1._rt));
                    X3(i,j,k) = c2 * (2._rt - tmp2 - 3._rt * tmp1 + 2._rt * tmp1 * tmp1 - 2._rt * I * om_mod * dt * tmp1)
                        / (4._rt * ep0 * om2_mod * (tmp1 - 1._rt));

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = tmp1_sqrt * (I - I * tmp2 * tmp2 - 2._rt * om_mod * dt) / (4._rt * ep0 * om_mod);
                }

                // Limits for nu = -omega_mod/omega
                if (nu == -om_mod/om) {

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = tmp2_sqrt * (1._rt - tmp1 * tmp1 + 2._rt * I * om_mod * dt) / (4._rt * ep0 * om2_mod);

                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * (- 3._rt + 4._rt * tmp1 - tmp1 * tmp1 - 2._rt * I * om_mod * dt)
                        / (4._rt * ep0 * om2_mod * (tmp1 - 1._rt));
                    X3(i,j,k) = c2 * (3._rt - 2._rt * tmp2 - 2._rt * tmp1 + tmp1 * tmp1 - 2._rt * I * om_mod * dt)
                        / (4._rt * ep0 * om2_mod * (tmp1 - 1._rt));

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = tmp2_sqrt * (- I + I * tmp1 * tmp1 - 2._rt * om_mod * dt) / (4._rt * ep0 * om_mod);
                }
            }

            // Limits for omega = 0 only
            else if (knorm_mod != 0. && knorm == 0.) {

                const amrex::Real om_mod  = c * knorm_mod;
                const amrex::Real om2_mod = om_mod * om_mod;

                C   (i,j,k) = std::cos(om_mod * dt);
                S_ck(i,j,k) = std::sin(om_mod * dt) / om_mod;
                T2(i,j,k) = 1._rt;

                // X1 multiplies i*(k \times J) in the update equation for B
                X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0 * om2_mod);

                // X2 multiplies rho_new in the update equation for E
                // X3 multiplies rho_old in the update equation for E
                X2(i,j,k) = c2 * (1._rt - S_ck(i,j,k) / dt) / (ep0 * om2_mod);
                X3(i,j,k) = c2 * (C(i,j,k) - S_ck(i,j,k) / dt) / (ep0 * om2_mod);

                // Coefficient multiplying J in update equation for E
                X4(i,j,k) = - S_ck(i,j,k) / ep0;

            }

            // Limits for omega_mod = 0 only
            else if (knorm_mod == 0. && knorm != 0.) {

                const amrex::Real om  = c * knorm;
                const amrex::Real om2 = om * om;
                const amrex::Real nu = - kv / om;
                const Complex theta      = amrex::exp(I * nu * om * dt * 0.5_rt);
                const Complex theta_star = amrex::exp(- I * nu * om * dt * 0.5_rt);

                C(i,j,k) = 1._rt;
                S_ck(i,j,k) = dt;
                T2(i,j,k) = theta * theta;

                if (nu != 0.) {
                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = (-theta_star + theta - I * nu * om * dt * theta)
                        / (ep0 * nu * nu * om2);

                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * (1._rt - T2(i,j,k) + I * nu * om * dt * T2(i,j,k)
                        + 0.5_rt * nu * nu * om2 * dt * dt * T2(i,j,k))
                        / (ep0 * nu * nu * om2 * (T2(i,j,k) - 1._rt));
                    X3(i,j,k) = c2 * (1._rt - T2(i,j,k) + I * nu * om * dt * T2(i,j,k)
                        + 0.5_rt * nu * nu * om2 * dt * dt)
                        / (ep0 * nu * nu * om2 * (T2(i,j,k) - 1._rt));

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = I * (theta - theta_star) / (ep0 * nu * om);
                }

                else {
                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = dt2 / (2._rt * ep0);

                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * dt2 / (6._rt * ep0);
                    X3(i,j,k) = - c2 * dt2 / (3._rt * ep0);

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = -dt / ep0;
                }
            }

            // Limits for omega_mod = 0 and omega = 0
            else if (knorm_mod == 0. && knorm == 0.) {

                C(i,j,k) = 1._rt;
                S_ck(i,j,k) = dt;
                T2(i,j,k) = 1._rt;

                // X1 multiplies i*(k \times J) in the update equation for B
                X1(i,j,k) = dt2 / (2._rt * ep0);

                // X2 multiplies rho_new in the update equation for E
                // X3 multiplies rho_old in the update equation for E
                X2(i,j,k) = c2 * dt2 / (6._rt * ep0);
                X3(i,j,k) = - c2 * dt2 / (3._rt * ep0);

                // Coefficient multiplying J in update equation for E
                X4(i,j,k) = -dt / ep0;
            }
        });
    }
}

void
ComovingPsatdAlgorithm::CurrentCorrection (const int lev,
                                           SpectralFieldData& field_data,
                                           std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                           const std::unique_ptr<amrex::MultiFab>& rho)
{
    // Profiling
    BL_PROFILE("ComovingAlgorithm::CurrentCorrection");

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
        const amrex::Real* const kx_arr = kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const amrex::Real* const modified_ky_arr = modified_ky_vec[mfi].dataPtr();
        const amrex::Real* const ky_arr = ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* const modified_kz_arr = modified_kz_vec[mfi].dataPtr();
        const amrex::Real* const kz_arr = kz_vec[mfi].dataPtr();

        // Local copy of member variables before GPU loop
        const amrex::Real dt = m_dt;

        // Comoving velocity
        const amrex::Real vx = m_v_comoving[0];
        const amrex::Real vy = m_v_comoving[1];
        const amrex::Real vz = m_v_comoving[2];

        // Loop over indices within one box
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Shortcuts for the values of J and rho
            const Complex Jx = fields(i,j,k,Idx::Jx);
            const Complex Jy = fields(i,j,k,Idx::Jy);
            const Complex Jz = fields(i,j,k,Idx::Jz);
            const Complex rho_old = fields(i,j,k,Idx::rho_old);
            const Complex rho_new = fields(i,j,k,Idx::rho_new);

            // k vector values, and coefficients
            const amrex::Real kx_mod = modified_kx_arr[i];
            const amrex::Real kx = kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const amrex::Real ky_mod = modified_ky_arr[j];
            const amrex::Real kz_mod = modified_kz_arr[k];
            const amrex::Real ky = ky_arr[j];
            const amrex::Real kz = kz_arr[k];
#else
            constexpr amrex::Real ky_mod = 0._rt;
            const     amrex::Real kz_mod = modified_kz_arr[j];
            constexpr amrex::Real ky = 0._rt;
            const     amrex::Real kz = kz_arr[j];
#endif
            constexpr Complex I = Complex{0._rt,1._rt};

            const amrex::Real knorm_mod = std::sqrt(kx_mod * kx_mod + ky_mod * ky_mod + kz_mod * kz_mod);

            // Correct J
            if (knorm_mod != 0._rt)
            {
                const Complex kmod_dot_J = kx_mod * Jx + ky_mod * Jy + kz_mod * Jz;
                const amrex::Real k_dot_v = kx * vx + ky * vy + kz * vz;

                if ( k_dot_v != 0._rt ) {

                    const Complex theta = amrex::exp(- I * k_dot_v * dt * 0.5_rt);
                    const Complex den = 1._rt - theta * theta;

                    fields(i,j,k,Idx::Jx) = Jx - (kmod_dot_J + k_dot_v * theta * (rho_new - rho_old) / den) * kx_mod / (knorm_mod * knorm_mod);
                    fields(i,j,k,Idx::Jy) = Jy - (kmod_dot_J + k_dot_v * theta * (rho_new - rho_old) / den) * ky_mod / (knorm_mod * knorm_mod);
                    fields(i,j,k,Idx::Jz) = Jz - (kmod_dot_J + k_dot_v * theta * (rho_new - rho_old) / den) * kz_mod / (knorm_mod * knorm_mod);

                } else {

                    fields(i,j,k,Idx::Jx) = Jx - (kmod_dot_J - I * (rho_new - rho_old) / dt) * kx_mod / (knorm_mod * knorm_mod);
                    fields(i,j,k,Idx::Jy) = Jy - (kmod_dot_J - I * (rho_new - rho_old) / dt) * ky_mod / (knorm_mod * knorm_mod);
                    fields(i,j,k,Idx::Jz) = Jz - (kmod_dot_J - I * (rho_new - rho_old) / dt) * kz_mod / (knorm_mod * knorm_mod);
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
ComovingPsatdAlgorithm::VayDeposition (const int /*lev*/,
                                       SpectralFieldData& /*field_data*/,
                                       std::array<std::unique_ptr<amrex::MultiFab>,3>& /*current*/)
{
    amrex::Abort("Vay deposition not implemented for comoving PSATD");
}

#endif // WARPX_USE_PSATD
