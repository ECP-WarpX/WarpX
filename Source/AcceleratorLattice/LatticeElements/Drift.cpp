/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Drift.H"

#include <AMReX_ParmParse.H>

#include <string>

Drift::Drift ()
    : LatticeElementBase("drift")
{
}

void
Drift::AddElement (amrex::ParmParse & pp_element, amrex::ParticleReal & z_location)
{
    AddElementBase(pp_element, z_location);
}
