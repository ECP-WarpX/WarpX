/* Copyright 2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_SPECIESPHYSICALPROPERTIES_H_
#define WARPX_SPECIESPHYSICALPROPERTIES_H_

#include "WarpXConst.H"
#include <AMReX_REAL.H>

#include <map>
#include <limits>

enum struct PhysicalSpecies{unspecified=0, electron, positron, photon, hydrogen};

std::map<std::string, PhysicalSpecies> SpeciesNames =
{
    { "unspecified", PhysicalSpecies::unspecified },
    { "electron", PhysicalSpecies::electron },
    { "positron", PhysicalSpecies::positron },
    { "photon",   PhysicalSpecies::photon },
    { "hydrogen", PhysicalSpecies::hydrogen }
};

std::map<PhysicalSpecies, amrex::Real> SpeciesMass =
{
    {PhysicalSpecies::unspecified, std::numeric_limits<amrex::Real>::quiet_NaN()},
    {PhysicalSpecies::electron, PhysConst::m_e},
    {PhysicalSpecies::positron, PhysConst::m_e},
    {PhysicalSpecies::photon,   std::numeric_limits<amrex::Real>::quiet_NaN()},
    {PhysicalSpecies::hydrogen, PhysConst::m_p}
};

std::map<PhysicalSpecies, amrex::Real> SpeciesCharge =
{
    {PhysicalSpecies::unspecified, std::numeric_limits<amrex::Real>::quiet_NaN()},
    {PhysicalSpecies::electron, -PhysConst::q_e},
    {PhysicalSpecies::positron,  PhysConst::q_e},
    {PhysicalSpecies::photon,    std::numeric_limits<amrex::Real>::quiet_NaN()},
    {PhysicalSpecies::hydrogen,  PhysConst::q_e}
};

#endif // WARPX_SPECIESPHYSICALPROPERTIES_H_
