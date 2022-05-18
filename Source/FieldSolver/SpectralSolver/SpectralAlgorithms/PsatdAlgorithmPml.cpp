/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PsatdAlgorithmPml.H"

#include "FieldSolver/SpectralSolver/SpectralFieldData.H"
#include "FieldSolver/SpectralSolver/SpectralKSpace.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpX_Complex.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_BaseFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_PODVector.H>

#include <cmath>

#if WARPX_USE_PSATD

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
PsatdAlgorithmPml::PsatdAlgorithmPml(const SpectralKSpace& spectral_kspace,
                                     const DistributionMapping& dm,
                                     const SpectralFieldIndex& spectral_index,
                                     const int norder_x, const int norder_y,
                                     const int norder_z, const bool nodal,
                                     const amrex::IntVect& fill_guards, const Real dt,
                                     const bool dive_cleaning, const bool divb_cleaning)
     // Initialize members of base class
     : SpectralBaseAlgorithm(spectral_kspace, dm, spectral_index, norder_x, norder_y, norder_z, nodal, fill_guards),
       m_spectral_index(spectral_index),
       m_dt(dt),
       m_dive_cleaning(dive_cleaning),
       m_divb_cleaning(divb_cleaning)
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
PsatdAlgorithmPml::pushSpectralFields(SpectralFieldData& f) const {

    const bool dive_cleaning = m_dive_cleaning;
    const bool divb_cleaning = m_divb_cleaning;

    const SpectralFieldIndex& Idx = m_spectral_index;

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
#if defined(WARPX_DIM_3D)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        const amrex::Real dt = m_dt;

        // Loop over indices within one box
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Record old values of the fields to be updated
            const Complex Exy = fields(i,j,k,Idx.Exy);
            const Complex Exz = fields(i,j,k,Idx.Exz);
            const Complex Eyx = fields(i,j,k,Idx.Eyx);
            const Complex Eyz = fields(i,j,k,Idx.Eyz);
            const Complex Ezx = fields(i,j,k,Idx.Ezx);
            const Complex Ezy = fields(i,j,k,Idx.Ezy);

            const Complex Bxy = fields(i,j,k,Idx.Bxy);
            const Complex Bxz = fields(i,j,k,Idx.Bxz);
            const Complex Byx = fields(i,j,k,Idx.Byx);
            const Complex Byz = fields(i,j,k,Idx.Byz);
            const Complex Bzx = fields(i,j,k,Idx.Bzx);
            const Complex Bzy = fields(i,j,k,Idx.Bzy);

            Complex Ex, Ey, Ez;
            Complex Bx, By, Bz;

            // Used only if dive_cleaning = true and divb_cleaning = true
            Complex Exx, Eyy, Ezz;
            Complex Bxx, Byy, Bzz;
            Complex Fx, Fy, Fz;
            Complex Gx, Gy, Gz;
            Complex F, G;

            if (!dive_cleaning && !divb_cleaning)
            {
                Ex = Exy + Exz;
                Ey = Eyx + Eyz;
                Ez = Ezx + Ezy;

                Bx = Bxy + Bxz;
                By = Byx + Byz;
                Bz = Bzx + Bzy;
            }
            else if (dive_cleaning && divb_cleaning)
            {
                Exx = fields(i,j,k,Idx.Exx);
                Eyy = fields(i,j,k,Idx.Eyy);
                Ezz = fields(i,j,k,Idx.Ezz);

                Bxx = fields(i,j,k,Idx.Bxx);
                Byy = fields(i,j,k,Idx.Byy);
                Bzz = fields(i,j,k,Idx.Bzz);

                Fx = fields(i,j,k,Idx.Fx);
                Fy = fields(i,j,k,Idx.Fy);
                Fz = fields(i,j,k,Idx.Fz);

                Gx = fields(i,j,k,Idx.Gx);
                Gy = fields(i,j,k,Idx.Gy);
                Gz = fields(i,j,k,Idx.Gz);

                Ex = Exx + Exy + Exz;
                Ey = Eyx + Eyy + Eyz;
                Ez = Ezx + Ezy + Ezz;

                Bx = Bxx + Bxy + Bxz;
                By = Byx + Byy + Byz;
                Bz = Bzx + Bzy + Bzz;

                F = Fx + Fy + Fz;
                G = Gx + Gy + Gz;
            }

            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];
#if defined(WARPX_DIM_3D)
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

                const amrex::Real C = C_arr(i,j,k);
                const amrex::Real S_ck = S_ck_arr(i,j,k);
                const amrex::Real inv_k2 = inv_k2_arr(i,j,k);

                const amrex::Real C1 = (kx2 * C + ky2 + kz2) * inv_k2;
                const amrex::Real C2 = (kx2 + ky2 * C + kz2) * inv_k2;
                const amrex::Real C3 = (kx2 + ky2 + kz2 * C) * inv_k2;
                const amrex::Real C4 = kx2 * (C - 1._rt) * inv_k2;
                const amrex::Real C5 = ky2 * (C - 1._rt) * inv_k2;
                const amrex::Real C6 = kz2 * (C - 1._rt) * inv_k2;
                const amrex::Real C7 = ky * kz * (1._rt - C) * inv_k2;
                const amrex::Real C8 = kx * kz * (1._rt - C) * inv_k2;
                const amrex::Real C9 = kx * ky * (1._rt - C) * inv_k2;

                if (!dive_cleaning && !divb_cleaning)
                {
                    const Complex C10 = I * c2 * kx * ky * kz * (dt - S_ck) * inv_k2;
                    const Complex C11 = I * c2 * ky2 * kz * (dt - S_ck) * inv_k2;
                    const Complex C12 = I * c2 * kz2 * ky * (dt - S_ck) * inv_k2;
                    const Complex C13 = I * c2 * kz2 * kx * (dt - S_ck) * inv_k2;
                    const Complex C14 = I * c2 * kx2 * kz * (dt - S_ck) * inv_k2;
                    const Complex C15 = I * c2 * kx2 * ky * (dt - S_ck) * inv_k2;
                    const Complex C16 = I * c2 * ky2 * kx * (dt - S_ck) * inv_k2;
                    const Complex C17 = I * c2 * kx * (ky2 * dt + (kz2 + kx2) * S_ck) * inv_k2;
                    const Complex C18 = I * c2 * kx * (kz2 * dt + (ky2 + kx2) * S_ck) * inv_k2;
                    const Complex C19 = I * c2 * ky * (kz2 * dt + (kx2 + ky2) * S_ck) * inv_k2;
                    const Complex C20 = I * c2 * ky * (kx2 * dt + (kz2 + ky2) * S_ck) * inv_k2;
                    const Complex C21 = I * c2 * kz * (kx2 * dt + (ky2 + kz2) * S_ck) * inv_k2;
                    const Complex C22 = I * c2 * kz * (ky2 * dt + (kx2 + kz2) * S_ck) * inv_k2;

                    const Complex C10_c2 = C10 / c2;
                    const Complex C11_c2 = C11 / c2;
                    const Complex C12_c2 = C12 / c2;
                    const Complex C13_c2 = C13 / c2;
                    const Complex C14_c2 = C14 / c2;
                    const Complex C15_c2 = C15 / c2;
                    const Complex C16_c2 = C16 / c2;
                    const Complex C17_c2 = C17 / c2;
                    const Complex C18_c2 = C18 / c2;
                    const Complex C19_c2 = C19 / c2;
                    const Complex C20_c2 = C20 / c2;
                    const Complex C21_c2 = C21 / c2;
                    const Complex C22_c2 = C22 / c2;

                    // Update E
                    fields(i,j,k,Idx.Exy) = C2  * Exy + C5  * Exz + C9  * Ey
                                          + C10 * Bx  + C11 * By  + C19 * Bz;

                    fields(i,j,k,Idx.Exz) = C6  * Exy + C3  * Exz + C8  * Ez
                                          - C10 * Bx  - C22 * By  - C12 * Bz;

                    fields(i,j,k,Idx.Eyz) = C3  * Eyz + C6  * Eyx + C7  * Ez
                                          + C21 * Bx  + C10 * By  + C13 * Bz;

                    fields(i,j,k,Idx.Eyx) = C9  * Ex + C4  * Eyz + C1  * Eyx
                                          - C14 * Bx - C10 * By  - C18 * Bz;

                    fields(i,j,k,Idx.Ezx) = C8  * Ex + C1  * Ezx + C4  * Ezy
                                          + C15 * Bx + C17 * By  + C10 * Bz;

                    fields(i,j,k,Idx.Ezy) = C7  * Ey + C5  * Ezx + C2  * Ezy
                                          - C20 * Bx - C16 * By  - C10 * Bz;

                    // Update B
                    fields(i,j,k,Idx.Bxy) = C2     * Bxy + C5     * Bxz + C9     * By
                                          - C10_c2 * Ex  - C11_c2 * Ey  - C19_c2 * Ez;

                    fields(i,j,k,Idx.Bxz) = C6     * Bxy + C3     * Bxz + C8     * Bz
                                          + C10_c2 * Ex  + C22_c2 * Ey  + C12_c2 * Ez;

                    fields(i,j,k,Idx.Byz) = C3     * Byz + C6     * Byx + C7     * Bz
                                          - C21_c2 * Ex  - C10_c2 * Ey  - C13_c2 * Ez;

                    fields(i,j,k,Idx.Byx) = C9     * Bx + C4     * Byz + C1     * Byx
                                          + C14_c2 * Ex + C10_c2 * Ey  + C18_c2 * Ez;

                    fields(i,j,k,Idx.Bzx) = C8     * Bx + C1     * Bzx + C4     * Bzy
                                          - C15_c2 * Ex - C17_c2 * Ey  - C10_c2 * Ez;

                    fields(i,j,k,Idx.Bzy) = C7     * By + C5     * Bzx + C2     * Bzy
                                          + C20_c2 * Ex + C16_c2 * Ey  + C10_c2 * Ez;
                }
                else if (dive_cleaning && divb_cleaning)
                {
                    const Complex C23 = I * c2 * kx * S_ck;
                    const Complex C24 = I * c2 * ky * S_ck;
                    const Complex C25 = I * c2 * kz * S_ck;

                    const Complex C23_c2 = C23 / c2;
                    const Complex C24_c2 = C24 / c2;
                    const Complex C25_c2 = C25 / c2;

                    // Update E
                    fields(i,j,k,Idx.Exx) = C1 * Exx + C4 * Exy + C4 * Exz
                                          - C9 * Ey  - C8 * Ez  + C23 * F;

                    fields(i,j,k,Idx.Exy) = C5 * Exx + C2  * Exy + C5 * Exz
                                          + C9 * Ey  + C24 * Bz  - C7 * G;

                    fields(i,j,k,Idx.Exz) = C6 * Exx + C6  * Exy + C3 * Exz
                                          + C8 * Ez  - C25 * By  + C7 * G;

                    fields(i,j,k,Idx.Eyx) = C9 * Ex  + C1  * Eyx + C4 * Eyy
                                          + C4 * Eyz - C23 * Bz  + C8 * G;

                    fields(i,j,k,Idx.Eyy) = - C9 * Ex  + C5 * Eyx + C2  * Eyy
                                            + C5 * Eyz - C7 * Ez  + C24 * F;

                    fields(i,j,k,Idx.Eyz) = C6 * Eyx + C6  * Eyy + C3 * Eyz
                                          + C7 * Ez  + C25 * Bx  - C8 * G;

                    fields(i,j,k,Idx.Ezx) = C8 * Ex  + C1  * Ezx + C4 * Ezy
                                          + C4 * Ezz + C23 * By  - C9 * G;

                    fields(i,j,k,Idx.Ezy) = C7 * Ey  + C5  * Ezx + C2 * Ezy
                                          + C5 * Ezz - C24 * Bx  + C9 * G;

                    fields(i,j,k,Idx.Ezz) = - C8 * Ex  - C7 * Ey  + C6  * Ezx
                                            + C6 * Ezy + C3 * Ezz + C25 * F;

                    // Update B
                    fields(i,j,k,Idx.Bxx) = C1 * Bxx + C4 * Bxy + C4     * Bxz
                                          - C9 * By  - C8 * Bz  + C23_c2 * G;

                    fields(i,j,k,Idx.Bxy) = - C24_c2 * Ez  + C5 * Bxx + C2 * Bxy
                                            + C5     * Bxz + C9 * By  + C7 * F;

                    fields(i,j,k,Idx.Bxz) = C25_c2 * Ey  + C6 * Bxx + C6 * Bxy
                                          + C3     * Bxz + C8 * Bz  - C7 * F;

                    fields(i,j,k,Idx.Byx) = C23_c2 * Ez  + C9 * Bx  + C1 * Byx
                                          + C4     * Byy + C4 * Byz - C8 * F;

                    fields(i,j,k,Idx.Byy) = - C9 * Bx  + C5 * Byx + C2     * Byy
                                            + C5 * Byz - C7 * Bz  + C24_c2 * G;

                    fields(i,j,k,Idx.Byz) = - C25_c2 * Ex  + C6 * Byx + C6 * Byy
                                            + C3     * Byz + C7 * Bz  + C8 * F;

                    fields(i,j,k,Idx.Bzx) = - C23_c2 * Ey  + C8 * Bx  + C1 * Bzx
                                            + C4     * Bzy + C4 * Bzz + C9 * F;

                    fields(i,j,k,Idx.Bzy) = C24_c2 * Ex  + C7 * By  + C5 * Bzx
                                          + C2     * Bzy + C5 * Bzz - C9 * F;

                    fields(i,j,k,Idx.Bzz) = - C8 * Bx  - C7 * By  + C6     * Bzx
                                            + C6 * Bzy + C3 * Bzz + C25_c2 * G;

                    // Update F
                    fields(i,j,k,Idx.Fx) = C23_c2 * Ex + C8 * By - C9 * Bz
                                         + C1     * Fx + C4 * Fy + C4 * Fz;

                    fields(i,j,k,Idx.Fy) = C24_c2 * Ey - C7 * Bx + C9 * Bz
                                         + C5     * Fx + C2 * Fy + C5 * Fz;

                    fields(i,j,k,Idx.Fz) = C25_c2 * Ez + C7 * Bx - C8 * By
                                         + C6     * Fx + C6 * Fy + C3 * Fz;

                    // Update G
                    fields(i,j,k,Idx.Gx) = - C8 * Ey + C9 * Ez + C23 * Bx
                                           + C1 * Gx + C4 * Gy + C4  * Gz;

                    fields(i,j,k,Idx.Gy) = C7 * Ex - C9 * Ez + C24 * By
                                         + C5 * Gx + C2 * Gy + C5  * Gz;

                    fields(i,j,k,Idx.Gz) = - C7 * Ex + C8 * Ey + C25 * Bz
                                           + C6 * Gx + C6 * Gy + C3  * Gz;
                }
            }
        });
    }
}

void PsatdAlgorithmPml::InitializeSpectralCoefficients (
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
#if defined(WARPX_DIM_3D)
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
#if defined(WARPX_DIM_3D)
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
PsatdAlgorithmPml::CurrentCorrection (SpectralFieldData& /*field_data*/)
{
    amrex::Abort("Current correction not implemented for PML PSATD");
}

void
PsatdAlgorithmPml::VayDeposition (SpectralFieldData& /*field_data*/)
{
    amrex::Abort("Vay deposition not implemented for PML PSATD");
}

#endif // WARPX_USE_PSATD
