/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "CollisionBase.H"

CollisionBase::CollisionBase (std::string collision_name)
{

    // read collision species
    amrex::ParmParse pp_collision_name(collision_name);
    pp_collision_name.getarr("species", m_species_names);

    // number of time steps between collisions
    m_ndt = 1;
    pp_collision_name.query("ndt", m_ndt);

}

