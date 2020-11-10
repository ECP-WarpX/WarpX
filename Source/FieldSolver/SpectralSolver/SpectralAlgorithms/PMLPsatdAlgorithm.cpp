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
     : SpectralBaseAlgorithm(spectral_kspace, dm, norder_x, norder_y, norder_z, nodal)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef    = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    C1_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C2_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C3_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C4_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C5_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C6_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C7_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C8_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C9_coef   = SpectralRealCoefficients(ba, dm, 1, 0);
    C10_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C11_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C12_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C13_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C14_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C15_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C16_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C17_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C18_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C19_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C20_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C21_coef  = SpectralRealCoefficients(ba, dm, 1, 0);
    C22_coef  = SpectralRealCoefficients(ba, dm, 1, 0);

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
        Array4<const Real> C1_arr  = C1_coef[mfi].array();
        Array4<const Real> C2_arr  = C2_coef[mfi].array();
        Array4<const Real> C3_arr  = C3_coef[mfi].array();
        Array4<const Real> C4_arr  = C4_coef[mfi].array();
        Array4<const Real> C5_arr  = C5_coef[mfi].array();
        Array4<const Real> C6_arr  = C6_coef[mfi].array();
        Array4<const Real> C7_arr  = C7_coef[mfi].array();
        Array4<const Real> C8_arr  = C8_coef[mfi].array();
        Array4<const Real> C9_arr  = C9_coef[mfi].array();
        Array4<const Real> C10_arr = C10_coef[mfi].array();
        Array4<const Real> C11_arr = C11_coef[mfi].array();
        Array4<const Real> C12_arr = C12_coef[mfi].array();
        Array4<const Real> C13_arr = C13_coef[mfi].array();
        Array4<const Real> C14_arr = C14_coef[mfi].array();
        Array4<const Real> C15_arr = C15_coef[mfi].array();
        Array4<const Real> C16_arr = C16_coef[mfi].array();
        Array4<const Real> C17_arr = C17_coef[mfi].array();
        Array4<const Real> C18_arr = C18_coef[mfi].array();
        Array4<const Real> C19_arr = C19_coef[mfi].array();
        Array4<const Real> C20_arr = C20_coef[mfi].array();
        Array4<const Real> C21_arr = C21_coef[mfi].array();
        Array4<const Real> C22_arr = C22_coef[mfi].array();

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
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif

            const Real k_norm = std::sqrt(kx*kx + ky*ky + kz*kz);

            constexpr Real c2 = PhysConst::c * PhysConst::c;

            const Complex I = Complex{0._rt, 1._rt};

            const Real C1  = C1_arr(i,j,k);
            const Real C2  = C2_arr(i,j,k);
            const Real C3  = C3_arr(i,j,k);
            const Real C4  = C4_arr(i,j,k);
            const Real C5  = C5_arr(i,j,k);
            const Real C6  = C6_arr(i,j,k);
            const Real C7  = C7_arr(i,j,k);
            const Real C8  = C8_arr(i,j,k);
            const Real C9  = C9_arr(i,j,k);
            const Real C10 = C10_arr(i,j,k);
            const Real C11 = C11_arr(i,j,k);
            const Real C12 = C12_arr(i,j,k);
            const Real C13 = C13_arr(i,j,k);
            const Real C14 = C14_arr(i,j,k);
            const Real C15 = C15_arr(i,j,k);
            const Real C16 = C16_arr(i,j,k);
            const Real C17 = C17_arr(i,j,k);
            const Real C18 = C18_arr(i,j,k);
            const Real C19 = C19_arr(i,j,k);
            const Real C20 = C20_arr(i,j,k);
            const Real C21 = C21_arr(i,j,k);
            const Real C22 = C22_arr(i,j,k);

            if (k_norm != 0.) {

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

    constexpr Real one = 1._rt;

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
        Array4<Real> C    = C_coef   [mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();
        Array4<Real> C1   = C1_coef  [mfi].array();
        Array4<Real> C2   = C2_coef  [mfi].array();
        Array4<Real> C3   = C3_coef  [mfi].array();
        Array4<Real> C4   = C4_coef  [mfi].array();
        Array4<Real> C5   = C5_coef  [mfi].array();
        Array4<Real> C6   = C6_coef  [mfi].array();
        Array4<Real> C7   = C7_coef  [mfi].array();
        Array4<Real> C8   = C8_coef  [mfi].array();
        Array4<Real> C9   = C9_coef  [mfi].array();
        Array4<Real> C10  = C10_coef [mfi].array();
        Array4<Real> C11  = C11_coef [mfi].array();
        Array4<Real> C12  = C12_coef [mfi].array();
        Array4<Real> C13  = C13_coef [mfi].array();
        Array4<Real> C14  = C14_coef [mfi].array();
        Array4<Real> C15  = C15_coef [mfi].array();
        Array4<Real> C16  = C16_coef [mfi].array();
        Array4<Real> C17  = C17_coef [mfi].array();
        Array4<Real> C18  = C18_coef [mfi].array();
        Array4<Real> C19  = C19_coef [mfi].array();
        Array4<Real> C20  = C20_coef [mfi].array();
        Array4<Real> C21  = C21_coef [mfi].array();
        Array4<Real> C22  = C22_coef [mfi].array();

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real kx = modified_kx[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky[j];
            const Real kz = modified_kz[k];
#else
            constexpr Real ky = 0.;
            const     Real kz = modified_kz[j];
#endif

            // Calculate norm of vector
            const Real k_norm = std::sqrt(kx*kx + ky*ky + kz*kz);
            const Real k2 = k_norm * k_norm;

            const Real kx2 = kx*kx;
            const Real ky2 = ky*ky;
            const Real kz2 = kz*kz;

            // Calculate coefficients
            constexpr Real c  = PhysConst::c;
            constexpr Real c2 = c*c;

            // Coefficients for k_norm = 0 do not need to be set
            if (k_norm != 0.) {

                const Real omega = c * k_norm;

                C(i,j,k)    = std::cos(omega * dt);
                S_ck(i,j,k) = std::sin(omega * dt) / omega;

                C1 (i,j,k) = (kx2 * C(i,j,k) + ky2 + kz2) / k2;
                C2 (i,j,k) = (kx2 + ky2 * C(i,j,k) + kz2) / k2;
                C3 (i,j,k) = (kx2 + ky2 + kz2 * C(i,j,k)) / k2;
                C4 (i,j,k) = kx2 * (C(i,j,k) - one) / k2;
                C5 (i,j,k) = ky2 * (C(i,j,k) - one) / k2;
                C6 (i,j,k) = kz2 * (C(i,j,k) - one) / k2;
                C7 (i,j,k) = ky * kz * (one - C(i,j,k)) / k2;
                C8 (i,j,k) = kx * kz * (one - C(i,j,k)) / k2;
                C9 (i,j,k) = kx * ky * (one - C(i,j,k)) / k2;
                C10(i,j,k) = c2 * kx * ky * kz * (dt - S_ck(i,j,k)) / k2;
                C11(i,j,k) = c2 * ky2 * kz * (dt - S_ck(i,j,k)) / k2;
                C12(i,j,k) = c2 * kz2 * ky * (dt - S_ck(i,j,k)) / k2;
                C13(i,j,k) = c2 * kz2 * kx * (dt - S_ck(i,j,k)) / k2;
                C14(i,j,k) = c2 * kx2 * kz * (dt - S_ck(i,j,k)) / k2;
                C15(i,j,k) = c2 * kx2 * ky * (dt - S_ck(i,j,k)) / k2;
                C16(i,j,k) = c2 * ky2 * kx * (dt - S_ck(i,j,k)) / k2;
                C17(i,j,k) = c2 * kx * (ky2 * dt + (kz2 + kx2) * S_ck(i,j,k)) / k2;
                C18(i,j,k) = c2 * kx * (kz2 * dt + (ky2 + kx2) * S_ck(i,j,k)) / k2;
                C19(i,j,k) = c2 * ky * (kz2 * dt + (kx2 + ky2) * S_ck(i,j,k)) / k2;
                C20(i,j,k) = c2 * ky * (kx2 * dt + (kz2 + ky2) * S_ck(i,j,k)) / k2;
                C21(i,j,k) = c2 * kz * (kx2 * dt + (ky2 + kz2) * S_ck(i,j,k)) / k2;
                C22(i,j,k) = c2 * kz * (ky2 * dt + (kx2 + kz2) * S_ck(i,j,k)) / k2;
            }
        });
    }
}

void
PMLPsatdAlgorithm::CurrentCorrection (SpectralFieldData& /*field_data*/,
                                      std::array<std::unique_ptr<amrex::MultiFab>,3>& /*current*/,
                                      const std::unique_ptr<amrex::MultiFab>& /*rho*/)
{
    amrex::Abort("Current correction not implemented for PML PSATD");
}

void
PMLPsatdAlgorithm::VayDeposition (SpectralFieldData& /*field_data*/,
                                  std::array<std::unique_ptr<amrex::MultiFab>,3>& /*current*/)
{
    amrex::Abort("Vay deposition not implemented for PML PSATD");
}

#endif // WARPX_USE_PSATD
