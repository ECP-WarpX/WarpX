/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

// Forward declaration needed
class MultiParticleContainer;

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
