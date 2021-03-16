/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PMLPsatdAlgorithm.H"
#include "Utils/WarpXConst.H"

#include <cmath>


#if WARPX_USE_PSATD

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
PMLPsatdAlgorithm::PMLPsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                                     const DistributionMapping& dm,
                                     const int norder_x, const int norder_y,
                                     const int norder_z, const bool nodal, const Real dt)
     // Initialize members of base class
     : SpectralBaseAlgorithm(spectral_kspace, dm, norder_x, norder_y, norder_z, nodal),
       m_dt(dt)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef      = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    inv_k2_coef = SpectralRealCoefficients(ba, dm, 1, 0);

    InitializeSpectralCoefficients(spectral_kspace, dm, dt);
}

/* Advance the E and B field in spectral space (stored in `f`)
 * over one time step */
void
PMLPsatdAlgorithm::pushSpectralFields(SpectralFieldData& f) const {

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        const Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();

        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
        Array4<const Real> inv_k2_arr = inv_k2_coef[mfi].array();

        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        const amrex::Real dt = m_dt;

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Record old values of the fields to be updated
            using Idx = SpectralPMLIndex;

            const Complex Exy = fields(i,j,k,Idx::Exy);
            const Complex Exz = fields(i,j,k,Idx::Exz);
            const Complex Eyx = fields(i,j,k,Idx::Eyx);
            const Complex Eyz = fields(i,j,k,Idx::Eyz);
            const Complex Ezx = fields(i,j,k,Idx::Ezx);
            const Complex Ezy = fields(i,j,k,Idx::Ezy);

            const Complex Bxy = fields(i,j,k,Idx::Bxy);
            const Complex Bxz = fields(i,j,k,Idx::Bxz);
            const Complex Byx = fields(i,j,k,Idx::Byx);
            const Complex Byz = fields(i,j,k,Idx::Byz);
            const Complex Bzx = fields(i,j,k,Idx::Bzx);
            const Complex Bzy = fields(i,j,k,Idx::Bzy);

            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0._rt;
            const Real kz = modified_kz_arr[j];
#endif
            constexpr Real c2 = PhysConst::c * PhysConst::c;

            const Complex I = Complex{0._rt, 1._rt};

            const Real kx2 = kx*kx;
            const Real ky2 = ky*ky;
            const Real kz2 = kz*kz;
            const Real k_norm = std::sqrt(kx2 + ky2 + kz2);

            if (k_norm != 0._rt) {

                const Real C = C_arr(i,j,k);
                const Real S_ck = S_ck_arr(i,j,k);
                const Real inv_k2 = inv_k2_arr(i,j,k);

                const Real C1  = (kx2 * C + ky2 + kz2) * inv_k2;
                const Real C2  = (kx2 + ky2 * C + kz2) * inv_k2;
                const Real C3  = (kx2 + ky2 + kz2 * C) * inv_k2;
                const Real C4  = kx2 * (C - 1._rt) * inv_k2;
                const Real C5  = ky2 * (C - 1._rt) * inv_k2;
                const Real C6  = kz2 * (C - 1._rt) * inv_k2;
                const Real C7  = ky * kz * (1._rt - C) * inv_k2;
                const Real C8  = kx * kz * (1._rt - C) * inv_k2;
                const Real C9  = kx * ky * (1._rt - C) * inv_k2;
                const Real C10 = c2 * kx * ky * kz * (dt - S_ck) * inv_k2;
                const Real C11 = c2 * ky2 * kz * (dt - S_ck) * inv_k2;
                const Real C12 = c2 * kz2 * ky * (dt - S_ck) * inv_k2;
                const Real C13 = c2 * kz2 * kx * (dt - S_ck) * inv_k2;
                const Real C14 = c2 * kx2 * kz * (dt - S_ck) * inv_k2;
                const Real C15 = c2 * kx2 * ky * (dt - S_ck) * inv_k2;
                const Real C16 = c2 * ky2 * kx * (dt - S_ck) * inv_k2;
                const Real C17 = c2 * kx * (ky2 * dt + (kz2 + kx2) * S_ck) * inv_k2;
                const Real C18 = c2 * kx * (kz2 * dt + (ky2 + kx2) * S_ck) * inv_k2;
                const Real C19 = c2 * ky * (kz2 * dt + (kx2 + ky2) * S_ck) * inv_k2;
                const Real C20 = c2 * ky * (kx2 * dt + (kz2 + ky2) * S_ck) * inv_k2;
                const Real C21 = c2 * kz * (kx2 * dt + (ky2 + kz2) * S_ck) * inv_k2;
                const Real C22 = c2 * kz * (ky2 * dt + (kx2 + kz2) * S_ck) * inv_k2;

                // Update E
                fields(i,j,k,Idx::Exy) = C2 * Exy + C5 * Exz + C9 * (Eyz + Eyx)
                                       + I*C10 * (Bxy + Bxz) + I*C11 * (Byz + Byx) + I*C19 * (Bzx + Bzy);

                fields(i,j,k,Idx::Exz) = C6 * Exy + C3 * Exz + C8 * (Ezx + Ezy)
                                       - I*C10 * (Bxy + Bxz) - I*C22 * (Byz + Byx) - I*C12 * (Bzx + Bzy);

                fields(i,j,k,Idx::Eyz) = C3 * Eyz + C6 * Eyx + C7 * (Ezx + Ezy)
                                       + I*C21 * (Bxy + Bxz) + I*C10 * (Byz + Byx) + I*C13 * (Bzx + Bzy);

                fields(i,j,k,Idx::Eyx) = C9 * (Exy + Exz) + C4 * Eyz + C1 * Eyx
                                       - I*C14 * (Bxy + Bxz) - I*C10 * (Byz + Byx) - I*C18 * (Bzx + Bzy);

                fields(i,j,k,Idx::Ezx) = C8 * (Exy + Exz) + C1 * Ezx + C4 * Ezy
                                       + I*C15 * (Bxy + Bxz) + I*C17 * (Byz + Byx) + I*C10 * (Bzx + Bzy);

                fields(i,j,k,Idx::Ezy) = C7 * (Eyz + Eyx) + C5 * Ezx + C2 * Ezy
                                       - I*C20 * (Bxy + Bxz) - I*C16 * (Byz + Byx) - I*C10 * (Bzx + Bzy);

                // Update B
                fields(i,j,k,Idx::Bxy) = C2 * Bxy + C5 * Bxz + C9 * (Byz + Byx)
                                       - I*C10/c2 * (Exy + Exz) - I*C11/c2 * (Eyz + Eyx) - I*C19/c2 * (Ezx + Ezy);

                fields(i,j,k,Idx::Bxz) = C6 * Bxy + C3 * Bxz + C8 * (Bzx + Bzy)
                                       + I*C10/c2 * (Exy + Exz) + I*C22/c2 * (Eyz + Eyx) + I*C12/c2 * (Ezx + Ezy);

                fields(i,j,k,Idx::Byz) = C3 * Byz + C6 * Byx + C7 * (Bzx + Bzy)
                                       - I*C21/c2 * (Exy + Exz) - I*C10/c2 * (Eyz + Eyx) - I*C13/c2 * (Ezx + Ezy);

                fields(i,j,k,Idx::Byx) = C9 * (Bxy + Bxz) + C4 * Byz + C1 * Byx
                                       + I*C14/c2 * (Exy + Exz) + I*C10/c2 * (Eyz + Eyx) + I*C18/c2 * (Ezx + Ezy);

                fields(i,j,k,Idx::Bzx) = C8 * (Bxy + Bxz) + C1 * Bzx + C4 * Bzy
                                       - I*C15/c2 * (Exy + Exz) - I*C17/c2 * (Eyz + Eyx) - I*C10/c2 * (Ezx + Ezy);

                fields(i,j,k,Idx::Bzy) = C7 * (Byz + Byx) + C5 * Bzx + C2 * Bzy
                                       + I*C20/c2 * (Exy + Exz) + I*C16/c2 * (Eyz + Eyx) + I*C10/c2 * (Ezx + Ezy);
            }
        });
    }
}

void PMLPsatdAlgorithm::InitializeSpectralCoefficients (
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const amrex::Real dt)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Fill them with the right values:
    // Loop over boxes and allocate the corresponding coefficients
    // for each box owned by the local MPI proc
    for (MFIter mfi(ba, dm); mfi.isValid(); ++mfi) {

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
        Array4<Real> inv_k2 = inv_k2_coef[mfi].array();

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real kx = modified_kx[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky[j];
            const Real kz = modified_kz[k];
#else
            constexpr Real ky = 0._rt;
            const     Real kz = modified_kz[j];
#endif

            // Calculate norm of vector
            const Real k_norm = std::sqrt(kx*kx + ky*ky + kz*kz);
            const Real k2 = k_norm * k_norm;

            // Calculate coefficients
            constexpr Real c  = PhysConst::c;

            // Coefficients for k_norm = 0 do not need to be set
            if (k_norm != 0._rt) {
                C(i,j,k) = std::cos(c * k_norm * dt);
                S_ck(i,j,k) = std::sin(c * k_norm * dt) / (c * k_norm);
                inv_k2(i,j,k) = 1._rt / k2;
            }
        });
    }
}

void
PMLPsatdAlgorithm::CurrentCorrection (const int /*lev*/,
                                      SpectralFieldData& /*field_data*/,
                                      std::array<std::unique_ptr<amrex::MultiFab>,3>& /*current*/,
                                      const std::unique_ptr<amrex::MultiFab>& /*rho*/)
{
    amrex::Abort("Current correction not implemented for PML PSATD");
}

void
PMLPsatdAlgorithm::VayDeposition (const int /*lev*/,
                                  SpectralFieldData& /*field_data*/,
                                  std::array<std::unique_ptr<amrex::MultiFab>,3>& /*current*/)
{
    amrex::Abort("Vay deposition not implemented for PML PSATD");
}

#endif // WARPX_USE_PSATD
