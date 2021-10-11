/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLE_UTILS_H_
#define WARPX_PARTICLE_UTILS_H_

#include "Particles/WarpXParticleContainer.H"

#include <AMReX_DenseBins.H>
#include <AMReX_Particles.H>

#include <AMReX_BaseFwd.H>

namespace ParticleUtils {

    /**
     * \brief Find the particles and count the particles that are in each cell. More specifically
     * this function returns an amrex::DenseBins object containing an offset array and a permutation
     * array which can be used to loop over all the cells in a tile and apply an algorithm to
     * particles of a given species present in each cell.
     * Note that this does *not* rearrange particle arrays.
     *
     * @param[in] lev the index of the refinement level.
     * @param[in] mfi the MultiFAB iterator.
     * @param[in] ptile the particle tile.
     */
    amrex::DenseBins<WarpXParticleContainer::ParticleType>
    findParticlesInEachCell( int const lev, amrex::MFIter const& mfi,
                             WarpXParticleContainer::ParticleTileType const& ptile);

    /**
     * \brief Generate random unit vector in 3 dimensions
     * https://mathworld.wolfram.com/SpherePointPicking.html
     *
     * @param[out] x x-component of resulting random vector
     * @param[out] y y-component of resulting random vector
     * @param[out] z z-component of resulting random vector
     * @param[in] engine the random-engine
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void getRandomVector ( amrex::Real& x, amrex::Real& y, amrex::Real& z,
                        amrex::RandomEngine const& engine )
    {
        using std::sqrt;
        using std::cos;
        using std::sin;
        using namespace amrex::literals;

        amrex::Real const theta = amrex::Random(engine) * 2.0_rt * MathConst::pi;
        z = 2.0_rt * amrex::Random(engine) - 1.0_rt;
        amrex::Real const xy = sqrt(1_rt - z*z);
        x = xy * cos(theta);
        y = xy * sin(theta);
    }


    /** \brief Function to perform scattering of a particle that results in a
     * random velocity vector with given magnitude. This is used in collision events.
     *
     * @param[in/out] ux, uy, uz colliding particle's velocity
     * @param[in] vp velocity magnitude of the colliding particle after collision.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void RandomizeVelocity ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
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
}

#endif // WARPX_PARTICLE_UTILS_H_
