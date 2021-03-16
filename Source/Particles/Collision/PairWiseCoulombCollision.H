/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_PAIRWISECOULOMBCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_PAIRWISECOULOMBCOLLISION_H_

#include "Particles/MultiParticleContainer.H"
#include "CollisionBase.H"
#include <AMReX_ParmParse.H>

class PairWiseCoulombCollision
    : public CollisionBase
{
public:

    PairWiseCoulombCollision (std::string collision_name);

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
    void doCoulombCollisionsWithinTile (
        int const lev, amrex::MFIter const& mfi,
        WarpXParticleContainer& species1,
        WarpXParticleContainer& species2);

private:

    amrex::Real m_CoulombLog;
    bool m_isSameSpecies;

};

#endif // WARPX_PARTICLES_COLLISION_PAIRWISECOULOMBCOLLISION_H_
