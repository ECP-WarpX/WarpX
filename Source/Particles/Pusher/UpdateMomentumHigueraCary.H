/* Copyright 2019 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_HIGUERACARY_H_
#define WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_HIGUERACARY_H_

#include "Utils/WarpXConst.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>


/**
 * \brief Push the particle's positions over one timestep,
 * given the value of its momenta `ux`, `uy`, `uz`
 */

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumHigueraCary(
    amrex::ParticleReal& ux, amrex::ParticleReal& uy, amrex::ParticleReal& uz,
    const amrex::ParticleReal Ex, const amrex::ParticleReal Ey, const amrex::ParticleReal Ez,
    const amrex::ParticleReal Bx, const amrex::ParticleReal By, const amrex::ParticleReal Bz,
    const amrex::ParticleReal q, const amrex::ParticleReal m, const amrex::Real dt )
{
    using namespace amrex::literals;

    // Constants
    const amrex::ParticleReal qmt = 0.5_prt*q*dt/m;
    constexpr amrex::ParticleReal invclight = 1._prt/PhysConst::c;
    constexpr amrex::ParticleReal invclightsq = 1._prt/(PhysConst::c*PhysConst::c);
    // Compute u_minus
    const amrex::ParticleReal umx = ux + qmt*Ex;
    const amrex::ParticleReal umy = uy + qmt*Ey;
    const amrex::ParticleReal umz = uz + qmt*Ez;
    // Compute gamma squared of u_minus
    amrex::ParticleReal gamma = 1._prt + (umx*umx + umy*umy + umz*umz)*invclightsq;
    // Compute beta and betam squared
    const amrex::ParticleReal betax = qmt*Bx;
    const amrex::ParticleReal betay = qmt*By;
    const amrex::ParticleReal betaz = qmt*Bz;
    const amrex::ParticleReal betam = betax*betax + betay*betay + betaz*betaz;
    // Compute sigma
    const amrex::ParticleReal sigma = gamma - betam;
    // Get u*
    const amrex::ParticleReal ust = (umx*betax + umy*betay + umz*betaz)*invclight;
    // Get new gamma inversed
    gamma = 1._prt/std::sqrt(0.5_prt*(sigma + std::sqrt(sigma*sigma + 4._prt*(betam + ust*ust)) ));
    // Compute t
    const amrex::ParticleReal tx = gamma*betax;
    const amrex::ParticleReal ty = gamma*betay;
    const amrex::ParticleReal tz = gamma*betaz;
    // Compute s
    const amrex::ParticleReal s = 1._prt/(1._prt+(tx*tx + ty*ty + tz*tz));
    // Compute um dot t
    const amrex::ParticleReal umt = umx*tx + umy*ty + umz*tz;
    // Compute u_plus
    const amrex::ParticleReal upx = s*( umx + umt*tx + umy*tz - umz*ty );
    const amrex::ParticleReal upy = s*( umy + umt*ty + umz*tx - umx*tz );
    const amrex::ParticleReal upz = s*( umz + umt*tz + umx*ty - umy*tx );
    // Get new u
    ux = upx + qmt*Ex + upy*tz - upz*ty;
    uy = upy + qmt*Ey + upz*tx - upx*tz;
    uz = upz + qmt*Ez + upx*ty - upy*tx;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_HIGUERACARY_H_
