/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef DISTANCETOEB_H_
#define DISTANCETOEB_H_

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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::RealVect interp_normal (int i, int j, int k, const amrex::Real W[AMREX_SPACEDIM][2],
                               amrex::Array4<const amrex::Real> const& phi,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) noexcept
{
#if (defined WARPX_DIM_3D)
    amrex::RealVect normal{0.0, 0.0, 0.0};

    normal[0] -= phi(i,   j  , k  ) * dxi[0] * W[1][0] * W[2][0];
    normal[0] += phi(i+1, j  , k  ) * dxi[0] * W[1][0] * W[2][0];
    normal[0] -= phi(i,   j+1, k  ) * dxi[0] * W[1][1] * W[2][0];
    normal[0] += phi(i+1, j+1, k  ) * dxi[0] * W[1][1] * W[2][0];
    normal[0] -= phi(i,   j  , k+1) * dxi[0] * W[1][0] * W[2][1];
    normal[0] += phi(i+1, j  , k+1) * dxi[0] * W[1][0] * W[2][1];
    normal[0] -= phi(i  , j+1, k+1) * dxi[0] * W[1][1] * W[2][1];
    normal[0] += phi(i+1, j+1, k+1) * dxi[0] * W[1][1] * W[2][1];

    normal[1] -= phi(i,   j  , k  ) * dxi[1] * W[0][0] * W[2][0];
    normal[1] += phi(i  , j+1, k  ) * dxi[1] * W[0][0] * W[2][0];
    normal[1] -= phi(i+1, j  , k  ) * dxi[1] * W[0][1] * W[2][0];
    normal[1] += phi(i+1, j+1, k  ) * dxi[1] * W[0][1] * W[2][0];
    normal[1] -= phi(i,   j  , k+1) * dxi[1] * W[0][0] * W[2][1];
    normal[1] += phi(i  , j+1, k+1) * dxi[1] * W[0][0] * W[2][1];
    normal[1] -= phi(i+1, j  , k+1) * dxi[1] * W[0][1] * W[2][1];
    normal[1] += phi(i+1, j+1, k+1) * dxi[1] * W[0][1] * W[2][1];

    normal[2] -= phi(i  , j  , k  ) * dxi[2] * W[0][0] * W[1][0];
    normal[2] += phi(i  , j  , k+1) * dxi[2] * W[0][0] * W[1][0];
    normal[2] -= phi(i+1, j  , k  ) * dxi[2] * W[0][1] * W[1][0];
    normal[2] += phi(i+1, j  , k+1) * dxi[2] * W[0][1] * W[1][0];
    normal[2] -= phi(i,   j+1, k  ) * dxi[2] * W[0][0] * W[1][1];
    normal[2] += phi(i  , j+1, k+1) * dxi[2] * W[0][0] * W[1][1];
    normal[2] -= phi(i+1, j+1, k  ) * dxi[2] * W[0][1] * W[1][1];
    normal[2] += phi(i+1, j+1, k+1) * dxi[2] * W[0][1] * W[1][1];
#elif (defined WARPX_DIM_XZ)
    amrex::RealVect normal{0.0, 0.0};

    normal[0] -= phi(i,   j  , k) * dxi[0] * W[1][0];
    normal[0] += phi(i+1, j  , k) * dxi[0] * W[1][0];
    normal[0] -= phi(i,   j+1, k) * dxi[0] * W[1][1];
    normal[0] += phi(i+1, j+1, k) * dxi[0] * W[1][1];

    normal[1] -= phi(i,   j  , k) * dxi[1] * W[0][0];
    normal[1] += phi(i  , j+1, k) * dxi[1] * W[0][0];
    normal[1] -= phi(i+1, j  , k) * dxi[1] * W[0][1];
    normal[1] += phi(i+1, j+1, k) * dxi[1] * W[0][1];
#else
    amrex::RealVect normal{0.0, 0.0};
    amrex::ignore_unused(i, j, k, W, phi, dxi);
    amrex::Abort("Error: interp_distance not yet implemented in RZ");
#endif
    return normal;
}
}

#endif // AMREX_USE_EB
#endif // DISTANCETOEB_H_
