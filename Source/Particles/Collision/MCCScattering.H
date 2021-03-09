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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void getRandomVector ( amrex::Real& x, amrex::Real& y, amrex::Real& z,
                       amrex::RandomEngine const& engine )
{
    // generate random unit vector in 3 dimensions
    // https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
    amrex::Real theta = amrex::Random(engine) * 2.0 * MathConst::pi;
    z = 2.0 * amrex::Random(engine) - 1.0;
    amrex::Real xy = std::sqrt(1 - z*z);
    x = xy * std::cos(theta);
    y = xy * std::sin(theta);
}

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
    amrex::Real x, y, z;
    // generate random unit vector for the new velocity direction
    getRandomVector(x, y, z, engine);

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
                      amrex::ParticleReal& uz,
                      const amrex::ParticleReal uCOM_x,
                      const amrex::ParticleReal uCOM_y,
                      const amrex::ParticleReal uCOM_z)
{
    // transform to COM frame, reverse particle velocity and transform back
    ux = -1.0 * ux + 2.0 * uCOM_x;
    uy = -1.0 * uy + 2.0 * uCOM_y;
    uz = -1.0 * uz + 2.0 * uCOM_z;
}


/** \brief Function to perform charge exchange of an ion with a neutral
 * particle.
 * @param[in/out] ux, uy, uz colliding particle's velocity
 * @param[in] ua_x, ua_y, ua_z velocity of the neutral particle.
 */

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ChargeExchange ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                      amrex::ParticleReal& uz,
                      const amrex::ParticleReal ua_x,
                      const amrex::ParticleReal ua_y,
                      const amrex::ParticleReal ua_z)
{
    // swap ion velocity for neutral velocity
    ux = ua_x;
    uy = ua_y;
    uz = ua_z;
}


/** \brief Function to perform scattering of an electron from an excitation
 * event. The target neutral is assumed stationary and the excited particle
 * is not tracked.
 * @param[in/out] ux, uy, uz colliding particle's velocity
 * @param[in] vp velocity magnitude of the colliding particle after collision.
 */

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void Excitation ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                  amrex::ParticleReal& uz,
                  const amrex::ParticleReal vp,
                  amrex::RandomEngine const& engine )
{
    amrex::Real x, y, z;
    // generate random unit vector for the new velocity direction
    getRandomVector(x, y, z, engine);

    // scale new vector to have the desired magnitude
    ux = x * vp;
    uy = y * vp;
    uz = z * vp;
}

#endif // WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_
