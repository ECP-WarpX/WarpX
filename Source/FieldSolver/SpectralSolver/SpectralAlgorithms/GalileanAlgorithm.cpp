#include "GalileanAlgorithm.H"
#include "Utils/WarpXConst.H"

#include <cmath>


#if WARPX_USE_PSATD

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
GalileanAlgorithm::GalileanAlgorithm(const SpectralKSpace& spectral_kspace,
                         const DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal,
                         const Array<Real, 3>& v_galilean,
                         const Real dt,
                         const bool update_with_rho)
     // Initialize members of base class
     : m_update_with_rho(update_with_rho),
       m_dt(dt),
       m_v_galilean(v_galilean),
       SpectralBaseAlgorithm(spectral_kspace, dm, norder_x, norder_y, norder_z, nodal)
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

    InitializeSpectralCoefficients(spectral_kspace, dm, v_galilean, dt);
};

/* Advance the E and B field in spectral space (stored in `f`)
 * over one time step */
void
GalileanAlgorithm::pushSpectralFields(SpectralFieldData& f) const{

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
        Array4<const Complex> Theta2_arr = Theta2_coef[mfi].array();

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
            constexpr Complex I = Complex{0,1};
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Complex X1 = X1_arr(i,j,k);
            const Complex X2 = X2_arr(i,j,k);
            const Complex X3 = X3_arr(i,j,k);
            const Complex X4 = X4_arr(i,j,k);
            const Complex T2 = Theta2_arr(i,j,k);

            // Update E (see the original Galilean article)

            if (update_with_rho) {

                fields(i,j,k,Idx::Ex) = T2*C*Ex_old
                            + T2*S_ck*c2*I*(ky*Bz_old - kz*By_old)
                            + X4*Jx - I*(X2*rho_new - T2*X3*rho_old)*kx;

                fields(i,j,k,Idx::Ey) = T2*C*Ey_old
                            + T2*S_ck*c2*I*(kz*Bx_old - kx*Bz_old)
                            + X4*Jy - I*(X2*rho_new - T2*X3*rho_old)*ky;

                fields(i,j,k,Idx::Ez) = T2*C*Ez_old
                            + T2*S_ck*c2*I*(kx*By_old - ky*Bx_old)
                            + X4*Jz - I*(X2*rho_new - T2*X3*rho_old)*kz;
            } else {

                Complex k_dot_J = kx * Jx + ky * Jy + kz * Jz;
                Complex k_dot_E = kx * Ex_old + ky * Ey_old + kz * Ez_old;

                fields(i,j,k,Idx::Ex) = T2 * C * Ex_old + I * T2 * S_ck * c2
                    * (ky * Bz_old - kz * By_old) + X4 * Jx + X2 * k_dot_E * kx + X3 * k_dot_J * kx;

                fields(i,j,k,Idx::Ey) = T2 * C * Ey_old + I * T2 * S_ck * c2
                    * (kz * Bx_old - kx * Bz_old) + X4 * Jy + X2 * k_dot_E * ky + X3 * k_dot_J * ky;

                fields(i,j,k,Idx::Ez) = T2 * C * Ez_old + I * T2 * S_ck * c2
                    * (kx * By_old - ky * Bx_old) + X4 * Jz + X2 * k_dot_E * kz + X3 * k_dot_J * kz;
            }

            // Update B (see the original Galilean article)
            // Note: here X1 is T2*x1/(ep0*c*c*k_norm*k_norm), where
            // x1 has the same definition as in the original paper

            fields(i,j,k,Idx::Bx) = T2*C*Bx_old
                        - T2*S_ck*I*(ky*Ez_old - kz*Ey_old)
                        +      X1*I*(ky*Jz     - kz*Jy);

            fields(i,j,k,Idx::By) = T2*C*By_old
                        - T2*S_ck*I*(kz*Ex_old - kx*Ez_old)
                        +      X1*I*(kz*Jx     - kx*Jz);

            fields(i,j,k,Idx::Bz) = T2*C*Bz_old
                        - T2*S_ck*I*(kx*Ey_old - ky*Ex_old)
                        +      X1*I*(kx*Jy     - ky*Jx);
        });
    }
};


void GalileanAlgorithm::InitializeSpectralCoefficients(const SpectralKSpace& spectral_kspace,
                                    const amrex::DistributionMapping& dm,
                                    const Array<Real, 3>& v_galilean,
                                    const amrex::Real dt)
{

    const bool update_with_rho = m_update_with_rho;

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
        Array4<Complex> X1 = X1_coef[mfi].array();
        Array4<Complex> X2 = X2_coef[mfi].array();
        Array4<Complex> X3 = X3_coef[mfi].array();
        Array4<Complex> X4 = X4_coef[mfi].array();
        Array4<Complex> Theta2 = Theta2_coef[mfi].array();
        // Extract reals (for portability on GPU)
        Real vx = v_galilean[0];
#if (AMREX_SPACEDIM==3)
        Real vy = v_galilean[1];
#endif
        Real vz = v_galilean[2];

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
            const Complex I{0.,1.};

            // Auxiliary coefficients
            Complex X2_old, X3_old;

            // Calculate dot product with galilean velocity
            const Real kv = modified_kx[i]*vx +
#if (AMREX_SPACEDIM==3)
                modified_ky[j]*vy + modified_kz[k]*vz;
#else
                modified_kz[j]*vz;
#endif

            if (k_norm != 0._rt) {

                C(i,j,k) = std::cos(c * k_norm * dt);
                S_ck(i,j,k) = std::sin(c * k_norm * dt) / (c * k_norm);

                const Real nu = kv/(k_norm*c);
                const Complex theta = amrex::exp( 0.5_rt*I*kv*dt );
                const Complex theta_star = amrex::exp( -0.5_rt*I*kv*dt );
                const Complex e_theta = amrex::exp( I*c*k_norm*dt );

                Theta2(i,j,k) = theta*theta;

                if ( (nu != 1.) && (nu != 0) ) {

                    // Note: the coefficients X1, X2, X3 do not correspond
                    // exactly to the original Galilean paper, but the
                    // update equation have been modified accordingly so that
                    // the expressions/ below (with the update equations)
                    // are mathematically equivalent to those of the paper.
                    Complex x1 = 1._rt/(1._rt-nu*nu) *
                        (theta_star - C(i,j,k)*theta + I*kv*S_ck(i,j,k)*theta);

                    // x1, above, is identical to the original paper
                    X1(i,j,k) = theta*x1/(ep0*c*c*k_norm*k_norm);

                    // The difference betwen X2 and X3 below, and those
                    // from the original paper is the factor ep0*k_norm*k_norm
                    if (update_with_rho) {

                        X2(i,j,k) = (x1 - theta * (1._rt - C(i,j,k))) / (theta_star - theta) / (ep0 * k_norm * k_norm);
                        X3(i,j,k) = (x1 - theta_star * (1._rt - C(i,j,k))) / (theta_star - theta) / (ep0 * k_norm * k_norm);

                    } else {

                        X2_old = (x1 - theta * (1._rt - C(i,j,k))) / (theta_star - theta);
                        X3_old = (x1 - theta_star * (1._rt - C(i,j,k))) / (theta_star - theta);

                        X2(i,j,k) = Theta2(i,j,k) * (X2_old - X3_old) / (k_norm * k_norm);
                        X3(i,j,k) = I * X2_old * (Theta2(i,j,k) - 1._rt) / (ep0 * k_norm * k_norm * kv);
                    }

                    X4(i,j,k) = I*kv*X1(i,j,k) - theta*theta*S_ck(i,j,k)/ep0;
                }

                if ( nu == 0) {

                    X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0*c*c*k_norm*k_norm);

                    if (update_with_rho) {

                        X2(i,j,k) = (1._rt - S_ck(i,j,k) / dt) / (ep0 * k_norm * k_norm);
                        X3(i,j,k) = (C(i,j,k) - S_ck(i,j,k) / dt) / (ep0 * k_norm * k_norm);

                    } else {

                        X2(i,j,k) = (1._rt - C(i,j,k)) / (k_norm * k_norm);
                        X3(i,j,k) = -dt * (1._rt - S_ck(i,j,k) / dt) / (ep0 * k_norm * k_norm);
                    }

                    X4(i,j,k) = -S_ck(i,j,k)/ep0;
                }

                if ( nu == 1.) {

                    X1(i,j,k) = (1._rt - e_theta*e_theta + 2._rt*I*c*k_norm*dt) / (4._rt*c*c*ep0*k_norm*k_norm);

                    if (update_with_rho) {

                        X2(i,j,k) = (3._rt - 4._rt * e_theta + e_theta * e_theta + 2._rt * I * c * k_norm * dt)
                            / (4._rt * ep0 * k_norm * k_norm * (1._rt - e_theta));

                        X3(i,j,k) = (3._rt - 2._rt / e_theta - 2._rt * e_theta + e_theta * e_theta
                            - 2._rt * I * c * k_norm * dt) / (4._rt * ep0 * (e_theta - 1._rt) * k_norm * k_norm);
                    }

                    else {

                        X2(i,j,k) = (1._rt - C(i,j,k)) * e_theta / (k_norm * k_norm);

                        X3(i,j,k) = (2._rt * dt * c * k_norm - I * e_theta * e_theta + 4._rt * I * e_theta - 3._rt * I)
                            / (4._rt * ep0 * c * k_norm * k_norm * k_norm);
                    }

                    X4(i,j,k) = I*(-1._rt + e_theta*e_theta + 2._rt*I*c*k_norm*dt) / (4._rt*ep0*c*k_norm);
                }
            // Handle k_norm = 0, by using the analytical limit
            } else {

                C(i,j,k) = 1._rt;
                S_ck(i,j,k) = dt;

                X1(i,j,k) = dt*dt/(2._rt * ep0);

                if (update_with_rho) {

                    X2(i,j,k) = c * c * dt * dt / (6._rt * ep0);
                    X3(i,j,k) = - c * c * dt * dt / (3._rt * ep0);

                } else {

                    X2(i,j,k) = c * c * dt *dt * 0.5_rt;
                    X3(i,j,k) = - c * c * dt * dt * dt / (6._rt * ep0);
                }

                X4(i,j,k) = -dt/ep0;

                Theta2(i,j,k) = 1._rt;
            }
        });
    }
}

void
GalileanAlgorithm::CurrentCorrection (SpectralFieldData& field_data,
                                      std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                      const std::unique_ptr<amrex::MultiFab>& rho) {
    // Profiling
    BL_PROFILE("GalileanAlgorithm::CurrentCorrection");

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
        const amrex::Real* const modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const amrex::Real* const modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const amrex::Real* const modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Local copy of member variables before GPU loop
        const amrex::Real dt = m_dt;

        // Galilean velocity
        const amrex::Real vgx = m_v_galilean[0];
        const amrex::Real vgy = m_v_galilean[1];
        const amrex::Real vgz = m_v_galilean[2];

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
            const amrex::Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const amrex::Real ky = modified_ky_arr[j];
            const amrex::Real kz = modified_kz_arr[k];
#else
            constexpr amrex::Real ky = 0._rt;
            const amrex::Real kz = modified_kz_arr[j];
#endif
            constexpr Complex I = Complex{0._rt,1._rt};

            const amrex::Real k_norm = std::sqrt(kx * kx + ky * ky + kz * kz);

            // Correct J
            if (k_norm != 0._rt)
            {
                const Complex k_dot_J = kx * Jx + ky * Jy + kz * Jz;
                const amrex::Real k_dot_vg = kx * vgx + ky * vgy + kz * vgz;

                if ( k_dot_vg != 0._rt ) {

                    const Complex rho_old_mod = rho_old * amrex::exp(I * k_dot_vg * dt);
                    const Complex den = 1._rt - amrex::exp(I * k_dot_vg * dt);

                    fields(i,j,k,Idx::Jx) = Jx - (k_dot_J - k_dot_vg * (rho_new - rho_old_mod) / den) * kx / (k_norm * k_norm);
                    fields(i,j,k,Idx::Jy) = Jy - (k_dot_J - k_dot_vg * (rho_new - rho_old_mod) / den) * ky / (k_norm * k_norm);
                    fields(i,j,k,Idx::Jz) = Jz - (k_dot_J - k_dot_vg * (rho_new - rho_old_mod) / den) * kz / (k_norm * k_norm);

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
#endif // WARPX_USE_PSATD
