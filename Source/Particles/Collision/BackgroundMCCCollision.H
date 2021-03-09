/* Copyright 2021 Roelof Groenewald
 *
 * This file is part of WarpX.
 *
 * License: ????
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_

#include "Particles/MultiParticleContainer.H"
#include "CollisionBase.H"
#include "CrossSectionHandler.H"
#include <AMReX_ParmParse.H>

class BackgroundMCCCollision
    : public CollisionBase
{
public:

    BackgroundMCCCollision (std::string collision_name);

    /** Perform the collisions
     *
     * @param cur_time Current time
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc) override;

    /** Perform all background collisions within a tile
     *
     * @param pti particle iterator
     * @param species1 pointer to species container
     *
     */
    void doBackgroundCollisionsWithinTile (
        WarpXParIter& pti,
        WarpXParticleContainer& species1);

private:

    amrex::Vector< CrossSectionHandler > m_scattering_processes;

    amrex::Real m_mass1;

    amrex::Real m_background_temperature;
    amrex::Real m_background_density;
    amrex::Real m_background_mass;
    amrex::Real m_total_collision_prob;
    amrex::Real m_nu_max;

};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_
