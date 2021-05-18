/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_

// _______________________________________________WarpX Forward declarations
#include "Particles/MultiParticleContainer.def"
// ______________________________________________

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include <string>

class CollisionBase
{
public:

    CollisionBase (std::string collision_name);

    virtual void doCollisions (amrex::Real /*cur_time*/, MultiParticleContainer* /*mypc*/ ){}

    CollisionBase(CollisionBase const &) = delete;
    CollisionBase(CollisionBase &&) = delete;
    CollisionBase & operator=(CollisionBase const &) = delete;

    virtual ~CollisionBase() = default;
protected:

    amrex::Vector<std::string> m_species_names;
    int m_ndt;

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_
