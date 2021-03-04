/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_COLLISIONHANDLER_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONHANDLER_H_

#include "CollisionBase.H"
#include <AMReX_Vector.H>

// Forward declaration needed
class MultiParticleContainer;

/* \brief CollisionHandler is a light weight class that contains the
 * list of collisions to be done.
 */
class CollisionHandler
{
public:
    CollisionHandler ();

    /* Perform all of the collisions */
    void doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc);

private:

    amrex::Vector<std::string> collision_names;
    amrex::Vector<std::string> collision_types;
    amrex::Vector< std::unique_ptr<CollisionBase> > allcollisions;

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONHANDLER_H_
