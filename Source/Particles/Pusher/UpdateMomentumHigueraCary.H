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
void UpdateMomentumHigueraCary (
    amrex::ParticleReal& ux, amrex::ParticleReal& uy, amrex::ParticleReal& uz,
    const amrex::ParticleReal Ex, const amrex::ParticleReal Ey, const amrex::ParticleReal Ez,
    const amrex::ParticleReal Bx, const amrex::ParticleReal By, const amrex::ParticleReal Bz,
    const amrex::Real q, const amrex::Real m, const amrex::Real dt )
{
    using namespace amrex;
    
    // Constants
    const Real qmt = 0.5*q*dt/m;
    constexpr Real invclight = 1./PhysConst::c;
    constexpr Real invclightsq = 1./(PhysConst::c*PhysConst::c);
    // Compute u_minus
    const Real umx = ux + qmt*Ex;
    const Real umy = uy + qmt*Ey;
    const Real umz = uz + qmt*Ez;
    // Compute gamma squared of u_minus
    Real gamma = 1. + (umx*umx + umy*umy + umz*umz)*invclightsq;
    // Compute beta and betam squared
    const Real betax = qmt*Bx;
    const Real betay = qmt*By;
    const Real betaz = qmt*Bz;
    const Real betam = betax*betax + betay*betay + betaz*betaz;
    // Compute sigma
    const Real sigma = gamma - betam;
    // Get u*
    const Real ust = (umx*betax + umy*betay + umz*betaz)*invclight;
    // Get new gamma inversed
    gamma = 1./std::sqrt(0.5*(sigma + std::sqrt(sigma*sigma + 4.*(betam + ust*ust)) ));
    // Compute t
    const Real tx = gamma*betax;
    const Real ty = gamma*betay;
    const Real tz = gamma*betaz;
    // Compute s
    const Real s = 1./(1.+(tx*tx + ty*ty + tz*tz));
    // Compute um dot t
    const Real umt = umx*tx + umy*ty + umz*tz;
    // Compute u_plus
    const Real upx = s*( umx + umt*tx + umy*tz - umz*ty );
    const Real upy = s*( umy + umt*ty + umz*tx - umx*tz );
    const Real upz = s*( umz + umt*tz + umx*ty - umy*tx );
    // Get new u
    ux = upx + qmt*Ex + upy*tz - upz*ty;
    uy = upy + qmt*Ey + upz*tx - upx*tz;
    uz = upz + qmt*Ez + upx*ty - upy*tx;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_HIGUERACARY_H_
