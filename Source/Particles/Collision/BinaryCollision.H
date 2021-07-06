/* Copyright 2020 David Grote, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_

#include "CollisionBase.H"
#include "Particles/WarpXParticleContainer.H"

#include "Particles/MultiParticleContainer_fwd.H"
#include "Particles/WarpXParticleContainer_fwd.H"

#include <AMReX_DenseBins.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <string>


/**
 * \brief A class from which specific binary collision functors acting within a single cell are
 *        derived.
 */
class BinaryCollisionFunc
{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;

public:
    /**
     * \brief Virtual operator() of the abstract BinaryCollisionFunc class
     */
    AMREX_GPU_HOST_DEVICE
    virtual void operator() (
        index_type const /*I1s*/, index_type const /*I1e*/,
        index_type const /*I2s*/, index_type const /*I2e*/,
        index_type* /*I1*/,       index_type* /*I2*/,
        amrex::Real* /*u1x*/, amrex::Real* /*u1y*/, amrex::Real* /*u1z*/,
        amrex::Real* /*u2x*/, amrex::Real* /*u2y*/, amrex::Real* /*u2z*/,
        amrex::Real const * /*w1*/, amrex::Real const * /*w2*/,
        amrex::Real const  /*q1*/, amrex::Real const  /*q2*/,
        amrex::Real const  /*m1*/, amrex::Real const  /*m2*/,
        amrex::Real const  /*dt*/, amrex::Real const /*dV*/,
        amrex::RandomEngine const& /*engine*/) const = 0;

    /**
     * \brief Virtual destructor of the abstract BinaryCollisionFunc class
     */
    virtual ~BinaryCollisionFunc () =  default;
};

class BinaryCollision
    : public CollisionBase
{
public:

    BinaryCollision (std::string collision_name);

    /** Perform the collisions
     *
     * @param lev AMR level of the tile
     * @param cur_time Current time
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc) override;

    /** Perform all binary collisions within a tile
     *
     * @param mfi iterator for multifab
     * @param species1/2 pointer to species container
     *
     */
    void doCollisionsWithinTile (
        int const lev, amrex::MFIter const& mfi,
        WarpXParticleContainer& species1,
        WarpXParticleContainer& species2);

private:

    bool m_isSameSpecies;
    std::unique_ptr<BinaryCollisionFunc> m_binary_collision_functor;

};

#endif // WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_
