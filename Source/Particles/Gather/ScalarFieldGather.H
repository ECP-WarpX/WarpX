/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef SCALARFIELDGATHER_H_
#define SCALARFIELDGATHER_H_

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
#else
    amrex::ignore_unused(xp, yp, zp, plo, dxi, i, j, k, W);
    amrex::Abort("Error: compute_weights not yet implemented in RZ");
#endif
}

/**
 * \brief Interpolate field value based on surrounding indices and weights.
 *
 * \param i, j, p                Indices of position on grid
 * \param W                      2D array of weights for each neighbouring node
 * \param scalar_field           Array4 of the scalar field, either full array or tile.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real interp_field (int i, int j, int k,
                               const amrex::Real W[AMREX_SPACEDIM][2],
                               amrex::Array4<const amrex::Real> const& scalar_field) noexcept
{
    amrex::Real value = 0;
#if (defined WARPX_DIM_3D)
    value += scalar_field(i,   j  , k  ) * W[0][0] * W[1][0] * W[2][0];
    value += scalar_field(i+1, j  , k  ) * W[0][1] * W[1][0] * W[2][0];
    value += scalar_field(i,   j+1, k  ) * W[0][0] * W[1][1] * W[2][0];
    value += scalar_field(i+1, j+1, k  ) * W[0][1] * W[1][1] * W[2][0];
    value += scalar_field(i,   j  , k+1) * W[0][0] * W[1][0] * W[2][1];
    value += scalar_field(i+1, j  , k+1) * W[0][1] * W[1][0] * W[2][1];
    value += scalar_field(i  , j+1, k+1) * W[0][0] * W[1][1] * W[2][1];
    value += scalar_field(i+1, j+1, k+1) * W[0][1] * W[1][1] * W[2][1];
#elif (defined WARPX_DIM_XZ)
    value += scalar_field(i,   j ,  k) * W[0][0] * W[1][0];
    value += scalar_field(i+1, j ,  k) * W[0][1] * W[1][0];
    value += scalar_field(i,   j+1, k) * W[0][0] * W[1][1];
    value += scalar_field(i+1, j+1, k) * W[0][1] * W[1][1];
#else
    amrex::ignore_unused(i, j, k, W, scalar_field);
    amrex::Abort("Error: interp_distance not yet implemented in RZ");
#endif
    return value;
}

/**
 * \brief Scalar field gather for a single particle.
 *
 * \param xp, yp, zp                Particle position coordinates
 * \param scalar_field              Array4 of the scalar field, either full array or tile.
 * \param dxi                       3D cell spacing
 * \param lo                        Index lower bounds of domain.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real doGatherScalarField (const amrex::ParticleReal xp,
                                 const amrex::ParticleReal yp,
                                 const amrex::ParticleReal zp,
                                 amrex::Array4<const amrex::Real> const& scalar_field,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& lo) noexcept
{
    // first find the weight of surrounding nodes to use during interpolation
    int ii, jj, kk;
    amrex::Real W[AMREX_SPACEDIM][2];
    compute_weights(xp, yp, zp, lo, dxi, ii, jj, kk, W);

    return interp_field(ii, jj, kk, W, scalar_field);
}
#endif // SCALARFIELDGATHER_H_
