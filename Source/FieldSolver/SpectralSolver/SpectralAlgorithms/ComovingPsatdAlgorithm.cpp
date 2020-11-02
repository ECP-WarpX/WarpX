#include "ComovingPsatdAlgorithm.H"
#include "Utils/WarpXConst.H"

#include <cmath>


#if WARPX_USE_PSATD

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
ComovingPsatdAlgorithm::ComovingPsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                                               const DistributionMapping& dm,
                                               const int norder_x, const int norder_y,
                                               const int norder_z, const bool nodal,
                                               const Array<Real, 3>& v_comoving,
                                               const Real dt,
                                               const bool update_with_rho)
     // Initialize members of base class
     : SpectralBaseAlgorithm(spectral_kspace, dm, norder_x, norder_y, norder_z, nodal),
       //kx_vec(spectral_kspace.getKComponent(dm,spectral_kspace.m_realspace_ba,0,false)),
       kx_vec(spectral_kspace.getModifiedKComponent(dm,0,-1,false)),
#if (AMREX_SPACEDIM==3)
       //ky_vec(spectral_kspace.getKComponent(dm,spectral_kspace.m_realspace_ba,1,false)),
       //kz_vec(spectral_kspace.getKComponent(dm,spectral_kspace.m_realspace_ba,2,false)),
       ky_vec(spectral_kspace.getModifiedKComponent(dm,1,-1,false)),
       kz_vec(spectral_kspace.getModifiedKComponent(dm,2,-1,false)),
#else
       //kz_vec(spectral_kspace.getKComponent(dm,spectral_kspace.m_realspace_ba,1,false)),
       kz_vec(spectral_kspace.getModifiedKComponent(dm,1,-1,false)),
#endif
       m_v_comoving(v_comoving),
       m_dt(dt),
       m_update_with_rho(update_with_rho)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    X1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X3_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X4_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    Theta2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

    InitializeSpectralCoefficients(spectral_kspace, dm, dt);
}

/* Advance the E and B field in spectral space (stored in `f`) over one time step */
void
ComovingPsatdAlgorithm::pushSpectralFields (SpectralFieldData& f) const
{
    const bool update_with_rho = m_update_with_rho;

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        const Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();

        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
        Array4<const Complex> X1_arr = X1_coef[mfi].array();
        Array4<const Complex> X2_arr = X2_coef[mfi].array();
        Array4<const Complex> X3_arr = X3_coef[mfi].array();
        Array4<const Complex> X4_arr = X4_coef[mfi].array();

        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0._rt;
            const     Real kz = modified_kz_arr[j];
#endif
            // Physical constant c**2 and imaginary unit
            constexpr Real c2   = PhysConst::c*PhysConst::c;
            constexpr Complex I = Complex{0._rt,1._rt};

            // The definition of these coefficients is explained in more detail
            // in the function InitializeSpectralCoefficients below
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Complex X1 = X1_arr(i,j,k);
            const Complex X2 = X2_arr(i,j,k);
            const Complex X3 = X3_arr(i,j,k);
            const Complex X4 = X4_arr(i,j,k);

            // The equations in the following are the update equations for B and E

            // Update E
            if (update_with_rho) {

                // Ex
                fields(i,j,k,Idx::Ex) = C*Ex_old
                            + S_ck*c2*I*(ky*Bz_old - kz*By_old)
                            + X4*Jx - I*(X2*rho_new - X3*rho_old)*kx;
                // Ey
                fields(i,j,k,Idx::Ey) = C*Ey_old
                            + S_ck*c2*I*(kz*Bx_old - kx*Bz_old)
                            + X4*Jy - I*(X2*rho_new - X3*rho_old)*ky;
                // Ez
                fields(i,j,k,Idx::Ez) = C*Ez_old
                            + S_ck*c2*I*(kx*By_old - ky*Bx_old)
                            + X4*Jz - I*(X2*rho_new - X3*rho_old)*kz;
            } else {
                // TODO Update E without rho
            }

            // Update B
            // Bx
            fields(i,j,k,Idx::Bx) = C*Bx_old - S_ck*I*(ky*Ez_old - kz*Ey_old) + X1*I*(ky*Jz - kz*Jy);
            // By
            fields(i,j,k,Idx::By) = C*By_old - S_ck*I*(kz*Ex_old - kx*Ez_old) + X1*I*(kz*Jx - kx*Jz);
            // Bz
            fields(i,j,k,Idx::Bz) = C*Bz_old - S_ck*I*(kx*Ey_old - ky*Ex_old) + X1*I*(kx*Jy - ky*Jx);
        });
    }
}

void ComovingPsatdAlgorithm::InitializeSpectralCoefficients (const SpectralKSpace& spectral_kspace,
                                                             const amrex::DistributionMapping& dm,
                                                             const amrex::Real dt)
{
    const bool update_with_rho = m_update_with_rho;

    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Loop over boxes and allocate the corresponding coefficients for each box
    for (MFIter mfi(ba, dm); mfi.isValid(); ++mfi) {

        const Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const Real* kx   = modified_kx_vec[mfi].dataPtr();
        const Real* kx_c = kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* ky   = modified_ky_vec[mfi].dataPtr();
        const Real* ky_c = ky_vec[mfi].dataPtr();
#endif
        const Real* kz   = modified_kz_vec[mfi].dataPtr();
        const Real* kz_c = kz_vec[mfi].dataPtr();

        // Extract arrays for the coefficients
        Array4<Real> C = C_coef[mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();
        Array4<Complex> X1 = X1_coef[mfi].array();
        Array4<Complex> X2 = X2_coef[mfi].array();
        Array4<Complex> X3 = X3_coef[mfi].array();
        Array4<Complex> X4 = X4_coef[mfi].array();
        Array4<Complex> T2 = Theta2_coef[mfi].array();

        // Extract Galilean velocity
        Real vx = m_v_comoving[0];
#if (AMREX_SPACEDIM==3)
        Real vy = m_v_comoving[1];
#endif
        Real vz = m_v_comoving[2];

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of vector
            const Real knorm = std::sqrt(
                std::pow(kx[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(ky[j], 2) +
                std::pow(kz[k], 2));
#else
                std::pow(kz[j], 2));
#endif
            // Calculate norm of vector
            const Real knorm_c = std::sqrt(
                std::pow(kx_c[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(ky_c[j], 2) +
                std::pow(kz_c[k], 2));
#else
                std::pow(kz_c[j], 2));
#endif
            // Physical constants c, c**2, and epsilon_0, and imaginary unit
            constexpr Real c    = PhysConst::c;
            constexpr Real c2   = c*c;
            constexpr Real ep0  = PhysConst::ep0;
            constexpr Complex I = Complex{0._rt,1._rt};

            // Auxiliary coefficients used when update_with_rho=false
            const Real dt2 = dt * dt;
            Complex X2_old, X3_old;

            // Calculate dot product of k vector with Galilean velocity
            const Real kv = kx_c[i]*vx +
#if (AMREX_SPACEDIM==3)
                ky_c[j]*vy + kz_c[k]*vz;
#else
                kz_c[j]*vz;
#endif

            if (knorm != 0. && knorm_c != 0.) {

                // Auxiliary coefficients
                const Real om  = c * knorm;
                const Real om2 = om * om;
                const Real om_c  = c * knorm_c;
                const Real om2_c = om_c * om_c;
                const Complex tmp1 = amrex::exp(  I * om * dt);
                const Complex tmp2 = amrex::exp(- I * om * dt);
                const Complex tmp1_sqrt = amrex::exp(  I * om * dt * 0.5);
                const Complex tmp2_sqrt = amrex::exp(- I * om * dt * 0.5);

                C   (i,j,k) = std::cos(om * dt);
                S_ck(i,j,k) = std::sin(om * dt) / om;

                const Real    nu = kv / om_c;
                const Complex theta      = amrex::exp(  I * nu * om_c * dt * 0.5_rt);
                const Complex theta_star = amrex::exp(- I * nu * om_c * dt * 0.5_rt);

                T2(i,j,k) = theta * theta;

                if ( (nu != om/om_c) && (nu != -om/om_c) && (nu != 0.) ) {

                    Complex x1 = om2_c / (om2 - nu * nu * om2_c)
                        * (theta_star - theta * C(i,j,k) + I * nu * om_c * theta * S_ck(i,j,k));

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = x1 / (ep0 * om2_c);

                    if (update_with_rho) {
                        // X2 multiplies rho_new in the update equation for E
                        // X3 multiplies rho_old in the update equation for E
                        X2(i,j,k) = c2 * (x1 * om2 - theta * (1._rt - C(i,j,k)) * om2_c)
                            / (theta_star - theta) / (ep0 * om2_c * om2);
                        X3(i,j,k) = c2 * (x1 * om2 - theta_star * (1._rt - C(i,j,k)) * om2_c)
                            / (theta_star - theta) / (ep0 * om2_c * om2);
                    } else {
                    }

                    // X4 multiplies J in the update equation for E
                    X4(i,j,k) = I * nu * om_c * X1(i,j,k) - theta * S_ck(i,j,k) / ep0;
                }

                // Limits for nu = 0
                if (nu == 0.) {

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0 * om2);

                    if (update_with_rho) {
                        // X2 multiplies rho_new in the update equation for E
                        // X3 multiplies rho_old in the update equation for E
                        X2(i,j,k) = c2 * (1._rt - S_ck(i,j,k) / dt) / (ep0 * om2);
                        X3(i,j,k) = c2 * (C(i,j,k) - S_ck(i,j,k) / dt) / (ep0 * om2);
                    } else {
                    }

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = - S_ck(i,j,k) / ep0;
                }

                // Limits for nu = omega/omega_c
                if (nu == om/om_c) {

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = tmp2_sqrt * (1._rt - tmp1 * tmp1 + 2._rt * I * om * dt) / (4._rt * ep0 * om2);

                    if (update_with_rho) {
                        // X2 multiplies rho_new in the update equation for E
                        // X3 multiplies rho_old in the update equation for E
                        X2(i,j,k) = c2 * (- 3._rt + 4._rt * tmp1 - tmp1 * tmp1 - 2._rt * I * om * dt)
                            / (4._rt * ep0 * om2 * (tmp1 - 1._rt));
                        X3(i,j,k) = c2 * (3._rt - 2._rt * tmp2 - 2._rt * tmp1 + tmp1 * tmp1 - 2._rt * I * om * dt)
                            / (4._rt * ep0 * om2 * (tmp1 - 1._rt));
                    } else {
                    }

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = tmp2_sqrt * (- I + I * tmp1 * tmp1 - 2._rt * om * dt) / (4._rt * ep0 * om);
                }

                // Limits for nu = -omega/omega_c
                if (nu == -om/om_c) {

                    // X1 multiplies i*(k \times J) in the update equation for B
                    X1(i,j,k) = tmp1_sqrt * (1._rt - tmp2 * tmp2 - 2._rt * I * om * dt) / (4._rt * ep0 * om2);

                    if (update_with_rho) {
                        // X2 multiplies rho_new in the update equation for E
                        // X3 multiplies rho_old in the update equation for E
                        X2(i,j,k) = c2 * (- 4._rt + 3._rt * tmp1 + tmp2 - 2._rt * I * om * dt * tmp1)
                            / (4._rt * ep0 * om2 * (tmp1 - 1._rt));
                        X3(i,j,k) = c2 * (2._rt - tmp2 - 3._rt * tmp1 + 2._rt * tmp1 * tmp1 - 2._rt * I * om * dt * tmp1)
                            / (4._rt * ep0 * om2 * (tmp1 - 1._rt));
                    } else {
                    }

                    // Coefficient multiplying J in update equation for E
                    X4(i,j,k) = tmp1_sqrt * (I - I * tmp2 * tmp2 - 2._rt * om * dt) / (4._rt * ep0 * om);
                }
            }

            // Limits for omega_c = 0 only
            else if (knorm != 0. && knorm_c == 0.) {

                const Real om  = c * knorm;
                const Real om2 = om * om;

                C   (i,j,k) = std::cos(om * dt);
                S_ck(i,j,k) = std::sin(om * dt) / om;
                T2(i,j,k) = 1._rt;

                // X1 multiplies i*(k \times J) in the update equation for B
                X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0 * om2);

                if (update_with_rho) {
                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * (1._rt - S_ck(i,j,k) / dt) / (ep0 * om2);
                    X3(i,j,k) = c2 * (C(i,j,k) - S_ck(i,j,k) / dt) / (ep0 * om2);
                } else {
                }

                // Coefficient multiplying J in update equation for E
                X4(i,j,k) = - S_ck(i,j,k) / ep0;

            }

            // Limits for omega = 0 only
            else if (knorm == 0. && knorm_c != 0.) {

                const Real om_c  = c * knorm_c;
                const Real om2_c = om_c * om_c;
                const Real nu = kv / om_c;
                const Complex theta      = amrex::exp(I * nu * om_c * dt * 0.5_rt);
                const Complex theta_star = amrex::exp(- I * nu * om_c * dt * 0.5_rt);

                C(i,j,k) = 1._rt;
                S_ck(i,j,k) = dt;
                T2(i,j,k) = theta * theta;

                // X1 multiplies i*(k \times J) in the update equation for B
                X1(i,j,k) = (-theta_star + theta - I * nu * om_c * dt * theta)
                    / (ep0 * nu * nu * om2_c);

                if (update_with_rho) {
                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * (1._rt - T2(i,j,k) + I * nu * om_c * dt * T2(i,j,k)
                        + 0.5_rt * nu * nu * om2_c * dt * dt * T2(i,j,k))
                        / (ep0 * nu * nu * om2_c * (T2(i,j,k) - 1._rt));
                    X3(i,j,k) = c2 * (1._rt - T2(i,j,k) + I * nu * om_c * dt * T2(i,j,k)
                        + 0.5_rt * nu * nu * om2_c * dt * dt)
                        / (ep0 * nu * nu * om2_c * (T2(i,j,k) - 1._rt));
                } else {
                }

                // Coefficient multiplying J in update equation for E
                X4(i,j,k) = I * (theta - theta_star) / (ep0 * nu * om_c);
            }

            // Limits for omega = 0 and omega_c = 0
            else if (knorm == 0. && knorm_c == 0.) {

                C(i,j,k) = 1._rt;
                S_ck(i,j,k) = dt;
                T2(i,j,k) = 1._rt;

                // X1 multiplies i*(k \times J) in the update equation for B
                X1(i,j,k) = dt2 / (2._rt * ep0);

                if (update_with_rho) {
                    // X2 multiplies rho_new in the update equation for E
                    // X3 multiplies rho_old in the update equation for E
                    X2(i,j,k) = c2 * dt2 / (6._rt * ep0);
                    X3(i,j,k) = - c2 * dt2 / (3._rt * ep0);
                } else {
                }

                // Coefficient multiplying J in update equation for E
                X4(i,j,k) = -dt / ep0;
            }
        });
    }
}

void
ComovingPsatdAlgorithm::CurrentCorrection (SpectralFieldData& field_data,
                                           std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                           const std::unique_ptr<amrex::MultiFab>& rho)
{
    // Profiling
    BL_PROFILE("ComovingAlgorithm::CurrentCorrection");

    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of J and rho
    field_data.ForwardTransform(*current[0], Idx::Jx, 0);
    field_data.ForwardTransform(*current[1], Idx::Jy, 0);
    field_data.ForwardTransform(*current[2], Idx::Jz, 0);
    field_data.ForwardTransform(*rho, Idx::rho_old, 0);
    field_data.ForwardTransform(*rho, Idx::rho_new, 1);

    // Loop over boxes
    for (amrex::MFIter mfi(field_data.fields); mfi.isValid(); ++mfi){

        const amrex::Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> fields = field_data.fields[mfi].array();

        // Extract pointers for the k vectors
        const amrex::Real* const modified_kx_arr   = modified_kx_vec[mfi].dataPtr();
        const amrex::Real* const kx_arr = kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const amrex::Real* const modified_ky_arr   = modified_ky_vec[mfi].dataPtr();
        const amrex::Real* const ky_arr = ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* const modified_kz_arr   = modified_kz_vec[mfi].dataPtr();
        const amrex::Real* const kz_arr = kz_vec[mfi].dataPtr();

        // Local copy of member variables before GPU loop
        const amrex::Real dt = m_dt;

        // Comoving velocity
        const amrex::Real vgx = m_v_comoving[0];
        const amrex::Real vgy = m_v_comoving[1];
        const amrex::Real vgz = m_v_comoving[2];

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
            const amrex::Real kx   = modified_kx_arr[i];
            const amrex::Real kx_c = kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const amrex::Real ky   = modified_ky_arr[j];
            const amrex::Real kz   = modified_kz_arr[k];
            const amrex::Real ky_c = ky_arr[j];
            const amrex::Real kz_c = kz_arr[k];
#else
            constexpr amrex::Real ky   = 0._rt;
            const     amrex::Real kz   = modified_kz_arr[j];
            constexpr amrex::Real ky_c = 0._rt;
            const     amrex::Real kz_c = kz_arr[j];
#endif
            constexpr Complex I = Complex{0._rt,1._rt};

            const amrex::Real k_norm = std::sqrt(kx * kx + ky * ky + kz * kz);

            // Correct J
            if (k_norm != 0._rt)
            {
                const Complex k_dot_J = kx * Jx + ky * Jy + kz * Jz;
                const amrex::Real k_dot_vg = kx_c * vgx + ky_c * vgy + kz_c * vgz;

                if ( k_dot_vg != 0._rt ) {

                    const Complex theta = amrex::exp(I * k_dot_vg * dt * 0.5);
                    const Complex den = 1._rt - theta * theta;

                    fields(i,j,k,Idx::Jx) = Jx - (k_dot_J - k_dot_vg * theta * (rho_new - rho_old) / den) * kx / (k_norm * k_norm);
                    fields(i,j,k,Idx::Jy) = Jy - (k_dot_J - k_dot_vg * theta * (rho_new - rho_old) / den) * ky / (k_norm * k_norm);
                    fields(i,j,k,Idx::Jz) = Jz - (k_dot_J - k_dot_vg * theta * (rho_new - rho_old) / den) * kz / (k_norm * k_norm);

                } else {

                    fields(i,j,k,Idx::Jx) = Jx - (k_dot_J - I * (rho_new - rho_old) / dt) * kx / (k_norm * k_norm);
                    fields(i,j,k,Idx::Jy) = Jy - (k_dot_J - I * (rho_new - rho_old) / dt) * ky / (k_norm * k_norm);
                    fields(i,j,k,Idx::Jz) = Jz - (k_dot_J - I * (rho_new - rho_old) / dt) * kz / (k_norm * k_norm);
                }
            }
        });
    }

    // Backward Fourier transform of J
    field_data.BackwardTransform(*current[0], Idx::Jx, 0);
    field_data.BackwardTransform(*current[1], Idx::Jy, 0);
    field_data.BackwardTransform(*current[2], Idx::Jz, 0);
}

void
ComovingPsatdAlgorithm::VayDeposition (SpectralFieldData& /*field_data*/,
                                  std::array<std::unique_ptr<amrex::MultiFab>,3>& /*current*/)
{
    amrex::Abort("Vay deposition not implemented for comoving PSATD");
}

#endif // WARPX_USE_PSATD
