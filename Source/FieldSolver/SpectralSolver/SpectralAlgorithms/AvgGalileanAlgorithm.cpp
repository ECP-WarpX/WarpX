#include "FieldSolver/SpectralSolver/SpectralAlgorithms/AvgGalileanAlgorithm.H"
#include "Utils/WarpXConst.H"
#include <cmath>

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
AvgGalileanAlgorithm::AvgGalileanAlgorithm(const SpectralKSpace& spectral_kspace,
                         const DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal,
                         const amrex::Array<amrex::Real,3>& v_galilean,
                         const Real dt)
     // Initialize members of base class
     : SpectralBaseAlgorithm( spectral_kspace, dm,
                              norder_x, norder_y, norder_z, nodal )
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);

    C1_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S1_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    C3_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S3_coef = SpectralRealCoefficients(ba, dm, 1, 0);

    Psi1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    Psi2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    Psi3_coef = SpectralComplexCoefficients(ba, dm, 1, 0);


    X1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X3_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    X4_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    Theta2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

    A1_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    A2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

    Rhoold_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    Rhonew_coef = SpectralComplexCoefficients(ba, dm, 1, 0);
    Jcoef_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

    InitializeSpectralCoefficients(spectral_kspace, dm, v_galilean, dt);

}

void AvgGalileanAlgorithm::InitializeSpectralCoefficients(
         const SpectralKSpace& spectral_kspace,
         const amrex::DistributionMapping& dm,
         const Array<Real, 3>& v_galilean,
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
        Array4<Real> C1 = C1_coef[mfi].array();
        Array4<Real> S1 = S1_coef[mfi].array();
        Array4<Real> C3 = C3_coef[mfi].array();
        Array4<Real> S3 = S3_coef[mfi].array();

        Array4<Complex> Psi1 = Psi1_coef[mfi].array();
        Array4<Complex> Psi2 = Psi2_coef[mfi].array();
        Array4<Complex> Psi3 = Psi3_coef[mfi].array();
        Array4<Complex> X1 = X1_coef[mfi].array();
        Array4<Complex> X2 = X2_coef[mfi].array();
        Array4<Complex> X3 = X3_coef[mfi].array();
        Array4<Complex> X4 = X4_coef[mfi].array();
        Array4<Complex> Theta2 = Theta2_coef[mfi].array();
        Array4<Complex> A1 = A1_coef[mfi].array();
        Array4<Complex> A2 = A2_coef[mfi].array();

        Array4<Complex> CRhoold = Rhoold_coef[mfi].array();
        Array4<Complex> CRhonew = Rhonew_coef[mfi].array();
        Array4<Complex> Jcoef   = Jcoef_coef[mfi].array();
        // Extract reals (for portability on GPU)

        Real vx = v_galilean[0];
        Real vy = v_galilean[1];
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
            constexpr Real c2 = PhysConst::c*PhysConst::c;
            constexpr Real ep0 = PhysConst::ep0;
            const Complex I{0.,1.};
            if (k_norm != 0){

                C(i,j,k) = std::cos(c*k_norm*dt);
                S_ck(i,j,k) = std::sin(c*k_norm*dt)/(c*k_norm);

                C1(i,j,k) = std::cos(0.5_rt*c*k_norm*dt);
                S1(i,j,k) = std::sin(0.5_rt*c*k_norm*dt);
                C3(i,j,k) = std::cos(1.5_rt*c*k_norm*dt);
                S3(i,j,k) = std::sin(1.5_rt*c*k_norm*dt);

                // Calculate dot product with galilean velocity
                const Real kv = modified_kx[i]*vx +
#if (AMREX_SPACEDIM==3)
                                 modified_ky[j]*vy +
                                 modified_kz[k]*vz;
#else
                                 modified_kz[j]*vz;
#endif

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

                    Complex C_rho = I* c2 /( (1._rt-theta*theta) * ep0);

                    Psi1(i,j,k) = theta * ((S1(i,j,k) + I*nu*C1(i,j,k))
                                  - Theta2(i,j,k) * (S3(i,j,k) + I*nu*C3(i,j,k))) /(c*k_norm*dt * (nu*nu - 1._rt));
                    Psi2(i,j,k) = theta * ((C1(i,j,k) - I*nu*S1(i,j,k))
                                  - Theta2(i,j,k) * (C3(i,j,k) - I*nu*S3(i,j,k))) /(c2*k_norm*k_norm*dt * (nu*nu - 1._rt));
                    Psi3(i,j,k) = I * theta * (1._rt - theta*theta) /(c*k_norm*dt*nu);

                    A1(i,j,k) = (Psi1(i,j,k)  - 1._rt + I * kv*Psi2(i,j,k)    )/ (c2* k_norm*k_norm * (nu*nu - 1._rt));
                    A2(i,j,k) = (Psi3(i,j,k) - Psi1(i,j,k)) / (c2*k_norm*k_norm);

                    CRhoold(i,j,k) = C_rho * (theta*theta * A1(i,j,k) - A2(i,j,k));
                    CRhonew(i,j,k) = C_rho * (A2(i,j,k) - A1(i,j,k));
                    Jcoef(i,j,k) = (I*kv*A1(i,j,k) + Psi2(i,j,k))/ep0;
                    // x1, above, is identical to the original paper
                    X1(i,j,k) = theta*x1/(ep0*c*c*k_norm*k_norm);
                    // The difference betwen X2 and X3 below, and those
                    // from the original paper is the factor ep0*k_norm*k_norm
                    X2(i,j,k) = (x1 - theta*(1._rt - C(i,j,k)))
                                /(theta_star-theta)/(ep0*k_norm*k_norm);
                    X3(i,j,k) = (x1 - theta_star*(1._rt - C(i,j,k)))
                                /(theta_star-theta)/(ep0*k_norm*k_norm);
                    X4(i,j,k) = I*kv*X1(i,j,k) - theta*theta*S_ck(i,j,k)/ep0;
                }
                if ( nu == 0) {
                    X1(i,j,k) = (1._rt - C(i,j,k)) / (ep0*c*c*k_norm*k_norm);
                    X2(i,j,k) = (1._rt - S_ck(i,j,k)/dt) / (ep0*k_norm*k_norm);
                    X3(i,j,k) = (C(i,j,k) - S_ck(i,j,k)/dt) / (ep0*k_norm*k_norm);
                    X4(i,j,k) = -S_ck(i,j,k)/ep0;

                    Psi1(i,j,k) = (-S1(i,j,k) + S3(i,j,k)) / (c*k_norm*dt);
                    Psi2(i,j,k) = (-C1(i,j,k) + C3(i,j,k)) / (c2*k_norm*k_norm*dt);
                    Psi3(i,j,k) = 1._rt;
                    A1(i,j,k) = (c*k_norm*dt + S1(i,j,k) - S3(i,j,k)) / (c*c2 * k_norm*k_norm*k_norm * dt);
                    A2(i,j,k) =  (c*k_norm*dt + S1(i,j,k) - S3(i,j,k)) / (c*c2 * k_norm*k_norm*k_norm * dt);
                    CRhoold(i,j,k) = 2._rt * I * S1(i,j,k)  * ( dt*C(i,j,k) - S_ck(i,j,k))
                                    / (c*k_norm*k_norm*k_norm*dt*dt*ep0);
                    CRhonew(i,j,k) =  - I * (c2* k_norm*k_norm * dt*dt - C1(i,j,k) + C3(i,j,k))
                                    / (c2 * k_norm*k_norm*k_norm*k_norm * ep0 * dt*dt);
                    Jcoef(i,j,k) = (-C1(i,j,k) + C3(i,j,k)) / (c2*ep0*k_norm*k_norm*dt);
                }
                if ( nu == 1.) {
                    X1(i,j,k) = (1._rt - e_theta*e_theta + 2._rt*I*c*k_norm*dt) / (4._rt*c*c*ep0*k_norm*k_norm);
                    X2(i,j,k) = (3._rt - 4._rt*e_theta + e_theta*e_theta + 2._rt*I*c*k_norm*dt) / (4._rt*ep0*k_norm*k_norm*(1._rt- e_theta));
                    X3(i,j,k) = (3._rt - 2._rt/e_theta - 2._rt*e_theta + e_theta*e_theta - 2._rt*I*c*k_norm*dt) / (4._rt*ep0*(e_theta - 1._rt)*k_norm*k_norm);
                    X4(i,j,k) = I*(-1._rt + e_theta*e_theta + 2._rt*I*c*k_norm*dt) / (4._rt*ep0*c*k_norm);
                }

            } else { // Handle k_norm = 0, by using the analytical limit
              C(i,j,k) = 1._rt;
              S_ck(i,j,k) = dt;
              C1(i,j,k) = 1._rt;
              S1(i,j,k) =  0._rt;
              C3(i,j,k) = 1._rt;
              S3(i,j,k) = 0._rt;

              X1(i,j,k) = dt*dt/(2._rt * ep0);
              X2(i,j,k) = c2*dt*dt/(6._rt * ep0);
              X3(i,j,k) = - c2*dt*dt/(3._rt * ep0);
              X4(i,j,k) = -dt/ep0;
              Theta2(i,j,k) = 1._rt;

              Psi1(i,j,k) = 1._rt;
              Psi2(i,j,k) = -dt;
              Psi3(i,j,k) = 1._rt;
              A1(i,j,k) = 13._rt * dt*dt /24._rt;
              A2(i,j,k) = 13._rt * dt*dt /24._rt;
              CRhoold(i,j,k) = -I*c2 * dt*dt / (3._rt * ep0);
              CRhonew(i,j,k) = -5._rt*I*c2 * dt*dt / (24._rt * ep0);
              Jcoef(i,j,k) = -dt/ep0;
            }

        });
    }
};

/* Advance the E and B field in spectral space (stored in `f`)
 * over one time step */
void
AvgGalileanAlgorithm::pushSpectralFields(SpectralFieldData& f) const{

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
        Array4<const Complex> Psi1_arr = Psi1_coef[mfi].array();
        Array4<const Complex> Psi2_arr = Psi2_coef[mfi].array();
        Array4<const Complex> Psi3_arr = Psi3_coef[mfi].array();
        Array4<const Real> C1_arr = C1_coef[mfi].array();
        Array4<const Real> S1_arr = S1_coef[mfi].array();
        Array4<const Real> C3_arr = C3_coef[mfi].array();
        Array4<const Real> S3_arr = S3_coef[mfi].array();


        Array4<const Complex> A1_arr = A1_coef[mfi].array();
        Array4<const Complex> A2_arr = A2_coef[mfi].array();
        Array4<const Complex> Rhonew_arr = Rhonew_coef[mfi].array();
        Array4<const Complex> Rhoold_arr = Rhoold_coef[mfi].array();
        Array4<const Complex> Jcoef_arr =Jcoef_coef[mfi].array();
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
            using Idx = SpectralAvgFieldIndex;

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

            const Complex Ex_avg = fields(i,j,k,Idx::Ex_avg);
            const Complex Ey_avg= fields(i,j,k,Idx::Ey_avg);
            const Complex Ez_avg = fields(i,j,k,Idx::Ez_avg);
            const Complex Bx_avg = fields(i,j,k,Idx::Bx_avg);
            const Complex By_avg = fields(i,j,k,Idx::By_avg);
            const Complex Bz_avg = fields(i,j,k,Idx::Bz_avg);
            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];

#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            constexpr Real c = PhysConst::c;
            constexpr Real c2 = PhysConst::c*PhysConst::c;
            constexpr Real inv_ep0 = 1._rt/PhysConst::ep0;
            constexpr Complex I = Complex{0,1};

            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);

            const Real C1 = C1_arr(i,j,k);
            const Real C3 = C3_arr(i,j,k);
            const Real S1 = S1_arr(i,j,k);
            const Real S3 = S3_arr(i,j,k);

            const Complex X1 = X1_arr(i,j,k);
            const Complex X2 = X2_arr(i,j,k);
            const Complex X3 = X3_arr(i,j,k);
            const Complex X4 = X4_arr(i,j,k);
            const Complex T2 = Theta2_arr(i,j,k);

            const Complex Psi1 = Psi1_arr(i,j,k);
            const Complex Psi2 = Psi2_arr(i,j,k);
            const Complex Psi3 = Psi3_arr(i,j,k);
            const Complex A1 = A1_arr(i,j,k);
            const Complex A2 = A2_arr(i,j,k);
            const Complex CRhoold= Rhoold_arr(i,j,k);
            const Complex CRhonew= Rhonew_arr(i,j,k);
            const Complex Jcoef = Jcoef_arr(i,j,k);


            //Update E (see the original Galilean article)
            fields(i,j,k,Idx::Ex) = T2*C*Ex_old
                        + T2*S_ck*c2*I*(ky*Bz_old - kz*By_old)
                        + X4*Jx - I*(X2*rho_new - T2*X3*rho_old)*kx;
            fields(i,j,k,Idx::Ey) = T2*C*Ey_old
                        + T2*S_ck*c2*I*(kz*Bx_old - kx*Bz_old)
                        + X4*Jy - I*(X2*rho_new - T2*X3*rho_old)*ky;
            fields(i,j,k,Idx::Ez) = T2*C*Ez_old
                        + T2*S_ck*c2*I*(kx*By_old - ky*Bx_old)
                        + X4*Jz - I*(X2*rho_new - T2*X3*rho_old)*kz;
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

//Update the averaged E,B fields in time on the interval [(n-1/2)dx, (n+1/2)dx]
            fields(i,j,k,Idx::Ex_avg) = Psi1*Ex_old
                          - Psi2*c2*I*(ky*Bz_old - kz*By_old)
                          + Jcoef*Jx + ( CRhonew * rho_new +  CRhoold*rho_old )*kx;
            fields(i,j,k,Idx::Ey_avg) =  Psi1*Ey_old
                          - Psi2*c2*I*(kz*Bx_old - kx*Bz_old)
                          + Jcoef*Jy +( CRhonew * rho_new +  CRhoold*rho_old )*ky;
            fields(i,j,k,Idx::Ez_avg) =  Psi1*Ez_old
                          - Psi2*c2*I*(kx*By_old - ky*Bx_old)
                          + Jcoef*Jz     + ( CRhonew * rho_new +  CRhoold*rho_old )*kz;

            fields(i,j,k,Idx::Bx_avg) =  Psi1*Bx_old
                          + I*Psi2*(ky*Ez_old - kz*Ey_old)
                          + A1*I*(ky*Jz     - kz*Jy)*inv_ep0;
            fields(i,j,k,Idx::By_avg) =  Psi1*By_old
                          + I*Psi2*(kz*Ex_old - kx*Ez_old)
                          + A1*I*(kz*Jx     - kx*Jz)*inv_ep0;
            fields(i,j,k,Idx::Bz_avg) =  Psi1*Bz_old
                          + I*Psi2*(kx*Ey_old - ky*Ex_old)
                          + A1*I*(kx*Jy     - ky*Jx)*inv_ep0;
                        });
    }
};
