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
 *    given the value of its momenta `ux`, `uy`, `uz` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdatePosition(amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
                    const amrex::ParticleReal ux, const amrex::ParticleReal uy, const amrex::ParticleReal uz,
                    const amrex::Real dt )
{
    using namespace amrex::literals;

    constexpr amrex::Real inv_c2 = 1._rt/(PhysConst::c*PhysConst::c);

    // Compute inverse Lorentz factor
    const amrex::Real inv_gamma = 1._rt/std::sqrt(1._rt + (ux*ux + uy*uy + uz*uz)*inv_c2);
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
