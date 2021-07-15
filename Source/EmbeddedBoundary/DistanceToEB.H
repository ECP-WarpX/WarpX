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
void compute_weights (const amrex::ParticleReal xp,
                      const amrex::ParticleReal yp,
                      const amrex::ParticleReal zp,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      int& i, int& j, int& k, amrex::Real W[AMREX_SPACEDIM][2]) noexcept
{
#if (defined WARPX_DIM_3D)
    amrex::Real x = (xp - plo[0]) * dxi[0];
    amrex::Real y = (yp - plo[1]) * dxi[1];
    amrex::Real z = (zp - plo[2]) * dxi[2];

    i = static_cast<int>(amrex::Math::floor(x));
    j = static_cast<int>(amrex::Math::floor(y));
    k = static_cast<int>(amrex::Math::floor(z));

    W[0][1] = x - i;
    W[1][1] = y - j;
    W[2][1] = z - k;

    W[0][0] = 1.0 - W[0][1];
    W[1][0] = 1.0 - W[1][1];
    W[2][0] = 1.0 - W[2][1];
#elif (defined WARPX_DIM_XZ)
    amrex::Real x = (xp - plo[0]) * dxi[0];
    amrex::Real z = (zp - plo[1]) * dxi[1];

    i = static_cast<int>(amrex::Math::floor(x));
    j = static_cast<int>(amrex::Math::floor(z));
    k = 0;

    W[0][1] = x - i;
    W[1][1] = z - j;

    W[0][0] = 1.0 - W[0][1];
    W[1][0] = 1.0 - W[1][1];

    amrex::ignore_unused(yp);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real interp_distance (int i, int j, int k, const amrex::Real W[AMREX_SPACEDIM][2],
                             amrex::Array4<const amrex::Real> const& phi) noexcept
{
    amrex::Real phi_value = 0;
#if (defined WARPX_DIM_3D)
    phi_value += phi(i,   j  , k  ) * W[0][0] * W[1][0] * W[2][0];
    phi_value += phi(i+1, j  , k  ) * W[0][1] * W[1][0] * W[2][0];
    phi_value += phi(i,   j+1, k  ) * W[0][0] * W[1][1] * W[2][0];
    phi_value += phi(i+1, j+1, k  ) * W[0][1] * W[1][1] * W[2][0];
    phi_value += phi(i,   j  , k+1) * W[0][0] * W[1][0] * W[2][1];
    phi_value += phi(i+1, j  , k+1) * W[0][1] * W[1][0] * W[2][1];
    phi_value += phi(i  , j+1, k+1) * W[0][0] * W[1][1] * W[2][1];
    phi_value += phi(i+1, j+1, k+1) * W[0][1] * W[1][1] * W[2][1];
#elif (defined WARPX_DIM_XZ)
    phi_value += phi(i,   j ,  k) * W[0][0] * W[1][0];
    phi_value += phi(i+1, j ,  k) * W[0][1] * W[1][0];
    phi_value += phi(i,   j+1, k) * W[0][0] * W[1][1];
    phi_value += phi(i+1, j+1, k) * W[0][1] * W[1][1];
#endif
    return phi_value;
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
#endif
    return normal;
}
}


#endif
