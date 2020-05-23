/* Copyright 2019 David Grote, Maxence Thevenet, Remi Lehe
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_VAY_H_
#define WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_VAY_H_

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>


/** \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumVay(
    amrex::ParticleReal& ux, amrex::ParticleReal& uy, amrex::ParticleReal& uz,
    const amrex::ParticleReal Ex, const amrex::ParticleReal Ey, const amrex::ParticleReal Ez,
    const amrex::ParticleReal Bx, const amrex::ParticleReal By, const amrex::ParticleReal Bz,
    const amrex::Real q, const amrex::Real m, const amrex::Real dt )
{
    // Constants
    const amrex::Real econst = q/m*dt;
    const amrex::Real bconst = 0.5*q/m*dt;
    constexpr amrex::Real invclight = 1./PhysConst::c;
    // Compute initial gamma
    const amrex::Real uxc = ux / PhysConst::c;
    const amrex::Real uyc = uy / PhysConst::c;
    const amrex::Real uzc = uz / PhysConst::c;
    const amrex::Real inv_gamma = 1./std::sqrt(1. + uxc*uxc + uyc*uyc + uzc*uzc);
    // Get tau
    const amrex::Real taux = bconst*Bx;
    const amrex::Real tauy = bconst*By;
    const amrex::Real tauz = bconst*Bz;
    const amrex::Real tausq = taux*taux+tauy*tauy+tauz*tauz;
    // Get U', gamma'^2
    const amrex::Real uxpr = ux + econst*Ex + (uy*tauz-uz*tauy)*inv_gamma;
    const amrex::Real uypr = uy + econst*Ey + (uz*taux-ux*tauz)*inv_gamma;
    const amrex::Real uzpr = uz + econst*Ez + (ux*tauy-uy*taux)*inv_gamma;
    const amrex::Real uxprc = uxpr / PhysConst::c;
    const amrex::Real uyprc = uypr / PhysConst::c;
    const amrex::Real uzprc = uzpr / PhysConst::c;
    const amrex::Real gprsq = (1. + uxprc*uxprc + uyprc*uyprc + uzprc*uzprc);
    // Get u*
    const amrex::Real ust = (uxpr*taux + uypr*tauy + uzpr*tauz)*invclight;
    // Get new gamma
    const amrex::Real sigma = gprsq-tausq;
    const amrex::Real gisq = 2./(sigma + std::sqrt(sigma*sigma + 4.*(tausq + ust*ust)) );
    // Get t, s
    const amrex::Real bg = bconst*std::sqrt(gisq);
    const amrex::Real tx = bg*Bx;
    const amrex::Real ty = bg*By;
    const amrex::Real tz = bg*Bz;
    const amrex::Real s = 1./(1.+tausq*gisq);
    // Get t.u'
    const amrex::Real tu = tx*uxpr + ty*uypr + tz*uzpr;
    // Get new U
    ux = s*(uxpr+tx*tu+uypr*tz-uzpr*ty);
    uy = s*(uypr+ty*tu+uzpr*tx-uxpr*tz);
    uz = s*(uzpr+tz*tu+uxpr*ty-uypr*tx);
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_VAY_H_
