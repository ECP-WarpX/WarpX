/* Copyright 2021 Roelof Groenewald
 *
 * This file is part of WarpX.
 *
 * License: ????
 */
#ifndef WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_
#define WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_

#include "Utils/WarpXConst.H"

#include <AMReX_Random.H>

/** \brief Function to perform elastic scattering of a particle in the lab
 * frame. The particle velocities transformed to the COM frame where a hard
 * sphere collision occurs. The resulting particle velocities are transformed
 * back to the lab frame and the input particle's velocity is updated.
 * @param[in/out] ux, uy, uz colliding particle's velocity
 * @param[in] uCOM_x, uCOM_y, uCOM_z velocity of the center of momentum frame.
 * @param[in] engine random number generator.
 */

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ElasticScattering ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                         amrex::ParticleReal& uz, amrex::ParticleReal& uCOM_x,
                         amrex::ParticleReal& uCOM_y, amrex::ParticleReal& uCOM_z,
                         amrex::RandomEngine const& engine )
{
    // generate random unit vector for the new velocity direction
    amrex::Real theta = amrex::Random() * 2.0 * MathConst::pi;
    amrex::Real z = 2.0 * amrex::Random() - 1.0;
    amrex::Real xy = std::sqrt(1 - z*z);
    amrex::Real x = xy * std::cos(theta);
    amrex::Real y = xy * std::sin(theta);

    // transform to center of momentum frame
    ux -= uCOM_x;
    uy -= uCOM_y;
    uz -= uCOM_z;

    // scale new vector to have the same magnitude as the original velocity
    amrex::Real mag = std::sqrt(ux*ux + uy*uy + uz*uz);
    ux = x * mag;
    uy = y * mag;
    uz = z * mag;

    // transform back to lab frame
    ux += uCOM_x;
    uy += uCOM_y;
    uz += uCOM_z;
}


/** \brief Function to perform back scattering of a particle in the lab
 * frame. The particle velocity is transformed to the COM frame where it is
 * reversed. The resulting particle velocities are then transformed back to the
 * lab frame and the input particle's velocity is updated.
 * @param[in/out] ux, uy, uz colliding particle's velocity
 * @param[in] uCOM_x, uCOM_y, uCOM_z velocity of the center of momentum frame.
 */

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void BackScattering ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                         amrex::ParticleReal& uz, amrex::ParticleReal& uCOM_x,
                         amrex::ParticleReal& uCOM_y, amrex::ParticleReal& uCOM_z)
{
    // transform to COM frame, reverse particle velocity and transform back
    ux = -1.0 * ux + 2.0 * uCOM_x;
    uy = -1.0 * uy + 2.0 * uCOM_y;
    uz = -1.0 * uz + 2.0 * uCOM_z;
}

#endif // WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_
