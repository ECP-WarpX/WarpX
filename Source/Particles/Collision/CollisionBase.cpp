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
    amrex::ParmParse pp(collision_name);
    pp.getarr("species", m_species_names);

    // number of time steps between collisions
    m_ndt = 1;
    pp.query("ndt", m_ndt);

}

