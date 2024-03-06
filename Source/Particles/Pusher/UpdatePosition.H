/* Copyright 2019 David Grote, Maxence Thevenet, Remi Lehe
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_
#define WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_

#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>



/** \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz`.
 *    This uses the standard leapfrog algorithm
 *    x^{n+1} - x^{n} = dt*u^{n+1/2}/gamma^{n+1/2}
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdatePosition(amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
                    const amrex::ParticleReal ux, const amrex::ParticleReal uy, const amrex::ParticleReal uz,
                    const amrex::Real dt )
{
    using namespace amrex::literals;

    constexpr amrex::ParticleReal inv_c2 = 1._prt/(PhysConst::c*PhysConst::c);

    // Compute inverse Lorentz factor
    const amrex::ParticleReal inv_gamma = 1._prt/std::sqrt(1._prt + (ux*ux + uy*uy + uz*uz)*inv_c2);
    // Update positions over one time step
#if (AMREX_SPACEDIM >= 2)
    x += ux * inv_gamma * dt;
#else
    amrex::ignore_unused(x);
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ) // RZ pushes particles in 3D
    y += uy * inv_gamma * dt;
#else
    amrex::ignore_unused(y);
#endif
    z += uz * inv_gamma * dt;
}

/** \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz`.
 *    The implicit version is the Crank-Nicolson scheme,
 *    x^{n+1} - x^{n} = dt*(u^{n+1} + u^{n})/(gamma^{n+1} + gamma^{n})
 *    See Eqs. 15 and 17 in Chen, JCP 407 (2020) 109228.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdatePositionImplicit(amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
                            const amrex::ParticleReal ux_n, const amrex::ParticleReal uy_n, const amrex::ParticleReal uz_n,
                            const amrex::ParticleReal ux, const amrex::ParticleReal uy, const amrex::ParticleReal uz,
                            const amrex::Real dt )
{
    using namespace amrex::literals;

    constexpr amrex::ParticleReal inv_c2 = 1._prt/(PhysConst::c*PhysConst::c);

    // Compute inverse Lorentz factor, the average of gamma at time levels n and n+1
    // The ux,uy,uz are the velocities at time level n+1/2
    const amrex::ParticleReal ux_np1 = 2._prt*ux - ux_n;
    const amrex::ParticleReal uy_np1 = 2._prt*uy - uy_n;
    const amrex::ParticleReal uz_np1 = 2._prt*uz - uz_n;
    const amrex::ParticleReal gamma_n = std::sqrt(1._prt + (ux_n*ux_n + uy_n*uy_n + uz_n*uz_n)*inv_c2);
    const amrex::ParticleReal gamma_np1 = std::sqrt(1._prt + (ux_np1*ux_np1 + uy_np1*uy_np1 + uz_np1*uz_np1)*inv_c2);
    const amrex::ParticleReal inv_gamma = 2.0_prt/(gamma_n + gamma_np1);

    // Update positions over one time step
#if (AMREX_SPACEDIM >= 2)
    x += ux * inv_gamma * dt;
#else
    amrex::ignore_unused(x);
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ) // RZ pushes particles in 3D
    y += uy * inv_gamma * dt;
#else
    amrex::ignore_unused(y);
#endif
    z += uz * inv_gamma * dt;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_
