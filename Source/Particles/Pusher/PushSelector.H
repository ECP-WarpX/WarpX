/* Copyright 2020 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_PUSHSELECTOR_H_
#define WARPX_PARTICLES_PUSHER_PUSHSELECTOR_H_

// Import low-level single-particle kernels
#include "Particles/Pusher/UpdateMomentumBoris.H"
#include "Particles/Pusher/UpdateMomentumBorisWithRadiationReaction.H"
#include "Particles/Pusher/UpdateMomentumHigueraCary.H"
#include "Particles/Pusher/UpdateMomentumVay.H"
#include "Particles/Pusher/UpdatePosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_REAL.H>

#include <limits>

/**
 * \brief Push position and momentum for a single particle
 *
 * \param GetPosition               A functor for returning the particle position.
 * \param SetPosition               A functor for setting the particle position.
 * \param copyAttribs               A functor for storing the old u and x
 * \param i                         The index of the particle to work on
 * \param ux, uy, uz                Particle momentum
 * \param Ex, Ey, Ez                Electric field on particles.
 * \param Bx, By, Bz                Magnetic field on particles.
 * \param ion_lev                   Ionization level of this particle (0 if ioniziation not on)
 * \param m                         Mass of this species.
 * \param q                         Charge of this species.
 * \param pusher_algo               0: Boris, 1: Vay, 2: HigueraCary
 * \param do_crr                    Whether to do the classical radiation reaction
 * \param do_copy                   Whether to copy the old x and u for the BTD
 * \param do_sync                   Whether to include quantum synchrotron radiation (QSR)
 * \param t_chi_max                 Cutoff chi for QSR
 * \param dt                        Time step size
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void doParticlePush(const GetParticlePosition& GetPosition,
                    const SetParticlePosition& SetPosition,
                    const CopyParticleAttribs& copyAttribs,
                    const long i,
                    amrex::ParticleReal& ux,
                    amrex::ParticleReal& uy,
                    amrex::ParticleReal& uz,
                    const amrex::ParticleReal Ex,
                    const amrex::ParticleReal Ey,
                    const amrex::ParticleReal Ez,
                    const amrex::ParticleReal Bx,
                    const amrex::ParticleReal By,
                    const amrex::ParticleReal Bz,
                    const int ion_lev,
                    const amrex::ParticleReal m,
                    const amrex::ParticleReal q,
                    const int pusher_algo,
                    const int do_crr,
                    const int do_copy,
#ifdef WARPX_QED
                    const int do_sync,
                    const amrex::Real t_chi_max,
#endif
                    const amrex::Real dt)
{
    if (do_copy) copyAttribs(i);
    if (do_crr) {
#ifdef WARPX_QED
        if (do_sync) {
            auto chi = QedUtils::chi_ele_pos(m*ux, m*uy, m*uz,
                                            Ex, Ey, Ez,
                                            Bx, By, Bz);
            if (chi < t_chi_max) {
                UpdateMomentumBorisWithRadiationReaction(ux, uy, uz,
                                                         Ex, Ey, Ez, Bx,
                                                         By, Bz, q, m, dt);
            }
            else {
                UpdateMomentumBoris( ux, uy, uz,
                                     Ex, Ey, Ez, Bx,
                                     By, Bz, q, m, dt);
            }
            amrex::ParticleReal x, y, z;
            GetPosition(i, x, y, z);
            UpdatePosition(x, y, z, ux, uy, uz, dt );
            SetPosition(i, x, y, z);
        } else {
            UpdateMomentumBorisWithRadiationReaction(ux, uy, uz,
                                                     Ex, Ey, Ez, Bx,
                                                     By, Bz, q, m, dt);
            amrex::ParticleReal x, y, z;
            GetPosition(i, x, y, z);
            UpdatePosition(x, y, z, ux, uy, uz, dt );
            SetPosition(i, x, y, z);
        }
#else
        amrex::ParticleReal qp = q;
        if (ion_lev) { qp *= ion_lev; }
        UpdateMomentumBorisWithRadiationReaction(ux, uy, uz,
                                                 Ex, Ey, Ez, Bx,
                                                 By, Bz, qp, m, dt);
        amrex::ParticleReal x, y, z;
        GetPosition(i, x, y, z);
        UpdatePosition(x, y, z, ux, uy, uz, dt );
        SetPosition(i, x, y, z);
#endif
    } else if (pusher_algo == ParticlePusherAlgo::Boris) {
        amrex::ParticleReal qp = q;
        if (ion_lev) { qp *= ion_lev; }
        UpdateMomentumBoris( ux, uy, uz,
                             Ex, Ey, Ez, Bx,
                             By, Bz, qp, m, dt);
        amrex::ParticleReal x, y, z;
        GetPosition(i, x, y, z);
        UpdatePosition(x, y, z, ux, uy, uz, dt );
        SetPosition(i, x, y, z);
    } else if (pusher_algo == ParticlePusherAlgo::Vay) {
        amrex::ParticleReal qp = q;
        if (ion_lev){ qp *= ion_lev; }
        UpdateMomentumVay( ux, uy, uz,
                           Ex, Ey, Ez, Bx,
                           By, Bz, qp, m, dt);
        amrex::ParticleReal x, y, z;
        GetPosition(i, x, y, z);
        UpdatePosition(x, y, z, ux, uy, uz, dt );
        SetPosition(i, x, y, z);
    } else if (pusher_algo == ParticlePusherAlgo::HigueraCary) {
        amrex::ParticleReal qp = q;
        if (ion_lev){ qp *= ion_lev; }
        UpdateMomentumHigueraCary( ux, uy, uz,
                                   Ex, Ey, Ez, Bx,
                                   By, Bz, qp, m, dt);
        amrex::ParticleReal x, y, z;
        GetPosition(i, x, y, z);
        UpdatePosition(x, y, z, ux, uy, uz, dt );
        SetPosition(i, x, y, z);
    } else {
        amrex::Abort("Unknown particle pusher");
    }
}

#endif // WARPX_PARTICLES_PUSHER_SELECTOR_H_
