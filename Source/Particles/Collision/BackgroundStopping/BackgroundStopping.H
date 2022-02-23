/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNSTOPPING_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNSTOPPING_H_

#include "Particles/MultiParticleContainer.H"
#include "Particles/Collision/CollisionBase.H"

#include <AMReX_REAL.H>

#include <string>

class BackgroundStopping final
    : public CollisionBase
{
public:
    BackgroundStopping (std::string collision_name);

    virtual ~BackgroundStopping () = default;

    /** Perform the stopping calculation
     *
     * @param cur_time Current time
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc) override;

    /** Perform the stopping calculation within a tile
     *
     * @param pti particle iterator
     *
     */
    void doBackgroundStoppingWithinTile (WarpXParIter& pti, amrex::Real dt,
                                         amrex::Real species_mass, amrex::Real species_charge);

private:

    amrex::Real m_background_temperature;
    amrex::Real m_background_density;
    amrex::Real m_background_mass;

};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNSTOPPING_H_
