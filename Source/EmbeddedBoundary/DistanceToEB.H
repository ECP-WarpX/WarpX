/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef DISTANCETOEB_H_
#define DISTANCETOEB_H_

#include "Utils/TextMsg.H"

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_Array.H>

#ifdef AMREX_USE_EB

namespace DistanceToEB
{

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dot_product (const amrex::RealVect& a, const amrex::RealVect& b) noexcept
{
    return AMREX_D_TERM(a[0]*b[0], + a[1]*b[1], + a[2]*b[2]);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void normalize (amrex::RealVect& a) noexcept
{
    amrex::Real inv_norm = 1.0/std::sqrt(dot_product(a,a));
    AMREX_D_DECL(a[0] *= inv_norm,
                 a[1] *= inv_norm,
                 a[2] *= inv_norm);
}



// This function calculates the normal vector using the nodal and cell-centered data.
// i,j,k are the index of the nearest node to the left of the point at which we interpolate.
// W are the interpolation weight for the left and right nodes (for the 0th component and 1st component respectively)
// ic,jc,kc are the index of the nearest cell-center to the left of the point at which we interpolate.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::RealVect interp_normal (int i, int j, int k, const amrex::Real W[AMREX_SPACEDIM][2],
                               int ic, int jc, int kc, const amrex::Real Wc[AMREX_SPACEDIM][2],
                               amrex::Array4<const amrex::Real> const& phi,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) noexcept
{

#if (defined WARPX_DIM_3D)
    amrex::RealVect normal{0.0, 0.0, 0.0};
    for (int iic = 0; iic < 2; ++iic) {
        for (int kk = 0; kk < 2; ++kk) {
            for (int jj=0; jj< 2; ++jj) {
                for (int ii = 0; ii < 2; ++ii) {
                    int icstart = ic + iic;
                    amrex::Real sign = (ii%2)*2. - 1.;
                    int wccomp = static_cast<int>(iic%2);
                    int w1comp = static_cast<int>(jj%2);
                    int w2comp = static_cast<int>(kk%2);
                    normal[0] += sign * phi(icstart + ii, j + jj, k + kk) * dxi[0] * Wc[0][wccomp] * W[1][w1comp] * W[2][w2comp];
                }
            }
        }
    }
    for (int iic = 0; iic < 2; ++iic) {
        for (int kk = 0; kk < 2; ++kk) {
            for (int ii=0; ii< 2; ++ii) {
                for (int jj = 0; jj < 2; ++jj) {
                    int jcstart = jc + iic;
                    amrex::Real sign = (jj%2)*2. - 1.;
                    int wccomp = static_cast<int>(iic%2);
                    int w1comp = static_cast<int>(ii%2);
                    int w2comp = static_cast<int>(kk%2);
                    normal[1] += sign * phi(i + ii, jcstart + jj, k + kk) * dxi[1] * W[0][w1comp] * Wc[1][wccomp] * W[2][w2comp];
                }
            }
        }
    }
    for (int iic = 0; iic < 2; ++iic) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int ii=0; ii< 2; ++ii) {
                for (int kk = 0; kk < 2; ++kk) {
                    int kcstart = kc + iic;
                    amrex::Real sign = (kk%2)*2. - 1.;
                    int wccomp = static_cast<int>(iic%2);
                    int w1comp = static_cast<int>(ii%2);
                    int w2comp = static_cast<int>(jj%2);
                    normal[2] += sign * phi(i + ii, j + jj, kcstart + kk) * dxi[2] * W[0][w1comp] * W[1][w2comp] * Wc[2][wccomp];
                }
            }
        }
    }

#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    amrex::RealVect normal{0.0, 0.0};
    for (int iic = 0; iic < 2; ++iic) {
        for (int jj=0; jj< 2; ++jj) {
            for (int ii = 0; ii < 2; ++ii) {
                int icstart = ic + iic;
                amrex::Real sign = (ii%2)*2. - 1.;
                int wccomp = static_cast<int>(iic%2);
                int w1comp = static_cast<int>(jj%2);
                normal[0] += sign * phi(icstart + ii, j + jj, k) * dxi[0] * Wc[0][wccomp] * W[1][w1comp];
            }
        }
    }
    for (int iic = 0; iic < 2; ++iic) {
        for (int ii=0; ii< 2; ++ii) {
            for (int jj = 0; jj < 2; ++jj) {
                int jcstart = jc + iic;
                amrex::Real sign = (jj%2)*2. - 1.;
                int wccomp = static_cast<int>(iic%2);
                int w1comp = static_cast<int>(ii%2);
                normal[1] += sign * phi(i + ii, jcstart + jj, k) * dxi[1] * W[0][w1comp] * Wc[1][wccomp];
            }
        }
    }
    amrex::ignore_unused(kc);

#else
    amrex::ignore_unused(i, j, k, ic, jc, kc, W, Wc, phi, dxi);
    amrex::RealVect normal{0.0, 0.0};
    WARPX_ABORT_WITH_MESSAGE("Error: interp_distance not yet implemented in 1D");

#endif
    return normal;
}
}

#endif // AMREX_USE_EB
#endif // DISTANCETOEB_H_
