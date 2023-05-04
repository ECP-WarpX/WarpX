/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PsatdAlgorithmGalileanPml.H"
#include "Utils/WarpXConst.H"

#include <cmath>

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
PsatdAlgorithmGalileanPml::PsatdAlgorithmGalileanPml (
    const SpectralKSpace& spectral_kspace,
    const DistributionMapping& dm,
    const SpectralFieldIndex& spectral_index,
    const int norder_x, const int norder_y,
    const int norder_z, const short grid_type,
    const amrex::Vector<amrex::Real>& v_galilean,
    const amrex::Real dt)
    // Initialize members of base class
    : SpectralBaseAlgorithm(spectral_kspace, dm, spectral_index, norder_x, norder_y, norder_z, grid_type)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    Theta2_coef = SpectralComplexCoefficients(ba, dm, 1, 0);

    InitializeSpectralCoefficients(spectral_kspace, dm, v_galilean, dt);
}

/* Advance the E and B field in spectral space (stored in `f`) over one time step */
void PsatdAlgorithmGalileanPml::pushSpectralFields (SpectralFieldData& f) const
{
    const SpectralFieldIndex& Idx = m_spectral_index;

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi)
    {
        const Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();
        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
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
            const Complex Ex_old = fields(i,j,k,Idx.Exy) \
                                 + fields(i,j,k,Idx.Exz);
            const Complex Ey_old = fields(i,j,k,Idx.Eyx) \
                                 + fields(i,j,k,Idx.Eyz);
            const Complex Ez_old = fields(i,j,k,Idx.Ezx) \
                                 + fields(i,j,k,Idx.Ezy);
            const Complex Bx_old = fields(i,j,k,Idx.Bxy) \
                                 + fields(i,j,k,Idx.Bxz);
            const Complex By_old = fields(i,j,k,Idx.Byx) \
                                 + fields(i,j,k,Idx.Byz);
            const Complex Bz_old = fields(i,j,k,Idx.Bzx) \
                                 + fields(i,j,k,Idx.Bzy);
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
            const Complex I = Complex{0,1};
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);
            const Complex T2 = Theta2_arr(i,j,k);

            // Update E
            fields(i,j,k,Idx.Exy) = T2*(C*fields(i,j,k,Idx.Exy) + S_ck*c2*I*ky*Bz_old);
            fields(i,j,k,Idx.Exz) = T2*(C*fields(i,j,k,Idx.Exz) - S_ck*c2*I*kz*By_old);
            fields(i,j,k,Idx.Eyz) = T2*(C*fields(i,j,k,Idx.Eyz) + S_ck*c2*I*kz*Bx_old);
            fields(i,j,k,Idx.Eyx) = T2*(C*fields(i,j,k,Idx.Eyx) - S_ck*c2*I*kx*Bz_old);
            fields(i,j,k,Idx.Ezx) = T2*(C*fields(i,j,k,Idx.Ezx) + S_ck*c2*I*kx*By_old);
            fields(i,j,k,Idx.Ezy) = T2*(C*fields(i,j,k,Idx.Ezy) - S_ck*c2*I*ky*Bx_old);
            // Update B
            fields(i,j,k,Idx.Bxy) = T2*(C*fields(i,j,k,Idx.Bxy) - S_ck*I*ky*Ez_old);
            fields(i,j,k,Idx.Bxz) = T2*(C*fields(i,j,k,Idx.Bxz) + S_ck*I*kz*Ey_old);
            fields(i,j,k,Idx.Byz) = T2*(C*fields(i,j,k,Idx.Byz) - S_ck*I*kz*Ex_old);
            fields(i,j,k,Idx.Byx) = T2*(C*fields(i,j,k,Idx.Byx) + S_ck*I*kx*Ez_old);
            fields(i,j,k,Idx.Bzx) = T2*(C*fields(i,j,k,Idx.Bzx) - S_ck*I*kx*Ey_old);
            fields(i,j,k,Idx.Bzy) = T2*(C*fields(i,j,k,Idx.Bzy) + S_ck*I*ky*Ex_old);
        });
    }
}

void PsatdAlgorithmGalileanPml::InitializeSpectralCoefficients (
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const amrex::Vector<amrex::Real>& v_galilean,
    const amrex::Real dt)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;
    // Fill them with the right values:
    // Loop over boxes and allocate the corresponding coefficients
    // for each box owned by the local MPI proc
    for (MFIter mfi(ba, dm); mfi.isValid(); ++mfi)
    {
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
        Array4<Complex> Theta2 = Theta2_coef[mfi].array();
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
            const Complex I{0.,1.};
            if (k_norm != 0)
            {
                C(i,j,k) = std::cos(c*k_norm*dt);
                S_ck(i,j,k) = std::sin(c*k_norm*dt)/(c*k_norm);
                // Calculate dot product with Galilean velocity
                const Real kv = modified_kx[i]*vx +
#if (AMREX_SPACEDIM==3)
                                 modified_ky[j]*vy +
                                 modified_kz[k]*vz;
#else
                                 modified_kz[j]*vz;
#endif
                Theta2(i,j,k) = amrex::exp( I*kv*dt );
            }
            else // Handle k_norm=0 by using the analytical limits
            {
                C(i,j,k) = 1.;
                S_ck(i,j,k) = dt;
                Theta2(i,j,k) = 1.;
            }
        });
    }
}

void PsatdAlgorithmGalileanPml::CurrentCorrection (SpectralFieldData& /*field_data*/)
{
    amrex::Abort("Current correction not implemented for PML PSATD");
}

void PsatdAlgorithmGalileanPml::VayDeposition (SpectralFieldData& /*field_data*/)
{
    amrex::Abort("Vay deposition not implemented for PML PSATD");
}
