/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLE_UTILS_H_
#define WARPX_PARTICLE_UTILS_H_

#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"

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
     * \brief Return (relativistic) particle energy given velocity and mass.
     *
     * @param[in] v2 square of particle speed ( v dot v)
     * @param[in] mass particle mass
     * @param[out] energy particle energy in eV
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void getEnergy ( amrex::ParticleReal const v2, amrex::ParticleReal const mass,
                     amrex::ParticleReal& gamma, amrex::ParticleReal& energy )
    {
        using std::sqrt;
        using namespace amrex::literals;

        // Compute c^2 since it is used repeatedly
        constexpr auto c2 = PhysConst::c * PhysConst::c;

        gamma = sqrt(1.0_prt + v2 / c2);
        energy = mass * c2 * (gamma - 1.0_prt) / PhysConst::q_e;
    }

    /**
     * \brief Return (relativistic) collision energy assuming the target (with
     * mass M) is stationary and the projectile is approaching with the
     * the given velocity and mass m.
     *
     * @param[in] v2 square of particle speed ( v dot v)
     * @param[in] m, M mass of projectile and target, respectively
     * @param[out] energy particle energy in eV
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void getCollisionEnergy ( amrex::ParticleReal const v2, amrex::ParticleReal const m,
                     amrex::ParticleReal const M, amrex::ParticleReal& gamma,
                     amrex::ParticleReal& energy )
    {
        using std::sqrt;
        using namespace amrex::literals;

        // Compute c^2 since it is used repeatedly
        constexpr auto c2 = PhysConst::c * PhysConst::c;

        const auto beta = v2 / c2;
        gamma = sqrt(1.0_prt + beta);
        const auto m_tot = M + gamma*m;
        energy = (
            sqrt(m_tot*m_tot - gamma*gamma*m*m*beta)*c2 - m*c2 - M*c2
        ) / PhysConst::q_e;
    }

    /**
     * \brief Return angles used to rotate the given particle velocity vector
     * such that v' = (0, 0, |v|).
     *
     * @param[in] vx, vy, vz velocity components
     * @param[out] theta, phi rotation angles to align velocity with the z-axis
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void getRotationAngles ( amrex::ParticleReal const vx, amrex::ParticleReal const vy,
                             amrex::ParticleReal const vz, amrex::ParticleReal& theta,
                             amrex::ParticleReal& phi )
    {
        using std::atan;
        using std::cos;
        using std::sin;

        // rotation to yz plane
        theta = atan(-vy / vx);
        // rotation to z axis
        phi = atan((vx * cos(theta) - vy * sin(theta)) / vz );
    }

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
    void getRandomVector ( amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
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
     * random velocity vector with given magnitude. This is used in isotropic
     * collision events.
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
