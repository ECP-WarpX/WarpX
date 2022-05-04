/* Copyright 2021 Luca Fedeli Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SpeciesPhysicalProperties.H"

#include "Utils/WarpXConst.H"

#include <AMReX_AmrCore.H>

#include <limits>
#include <map>
#include <string>

using namespace species;

namespace {
    struct Properties
    {
        amrex::Real mass;
        amrex::Real charge;
    };

    const auto string_to_species = std::map<std::string, PhysicalSpecies>{
        {"unspecified", PhysicalSpecies::unspecified},
        {"electron"   , PhysicalSpecies::electron},
        {"positron"   , PhysicalSpecies::positron},
        {"photon"     , PhysicalSpecies::photon},
        {"hydrogen"   , PhysicalSpecies::hydrogen},
        {"proton"     , PhysicalSpecies::hydrogen},
        {"helium"     , PhysicalSpecies::helium},
        {"alpha"      , PhysicalSpecies::helium},
        {"boron"      , PhysicalSpecies::boron},
        {"boron10"    , PhysicalSpecies::boron10},
        {"boron11"    , PhysicalSpecies::boron11},
        {"carbon"     , PhysicalSpecies::carbon},
        {"nitrogen"   , PhysicalSpecies::nitrogen},
        {"oxygen"     , PhysicalSpecies::oxygen},
        {"argon"      , PhysicalSpecies::argon},
        {"copper"     , PhysicalSpecies::copper},
        {"xenon"      , PhysicalSpecies::xenon}
    };

    const auto species_to_string = std::map<PhysicalSpecies, std::string>{
        {PhysicalSpecies::unspecified, "unspecified"},
        {hysicalSpecies::electron    , "electron"},
        {PhysicalSpecies::positron   , "positron"},
        {PhysicalSpecies::photon     , "photon"},
        {PhysicalSpecies::hydrogen   , "hydrogen"},
        {PhysicalSpecies::helium     , "helium"},
        {PhysicalSpecies::boron      , "boron"},
        {PhysicalSpecies::boron10    , "boron10"},
        {PhysicalSpecies::boron11    , "boron11"},
        {PhysicalSpecies::carbon     , "carbon"},
        {PhysicalSpecies::nitrogen   , "nitrogen"},
        {PhysicalSpecies::oxygen     , "oxygen"},
        {PhysicalSpecies::argon      , "argon"},
        {PhysicalSpecies::copper     , "copper"},
        {PhysicalSpecies::xenon      , "xenon"}
    };

    constexpr auto quiet_NaN = std::numeric_limits<amrex::Real>::quiet_NaN();

    const
    std::map<PhysicalSpecies,Properties> species_to_properties
    {
        {PhysicalSpecies::unspecified,
            Properties{quiet_NaN, quiet_NaN}},
        {PhysicalSpecies::electron,
            Properties{PhysConst::m_e, -PhysConst::q_e}},
        {PhysicalSpecies::positron,
            Properties{PhysConst::m_e, PhysConst::q_e}},
        {PhysicalSpecies::photon ,
            Properties{amrex::Real(0.0), PhysConst::q_e}},
        {PhysicalSpecies::hydrogen,
            Properties{PhysConst::m_p, PhysConst::q_e}},
        {PhysicalSpecies::helium,
            Properties{PhysConst::m_p * amrex::Real(3.97369), PhysConst::q_e * amrex::Real(2.0)}},
        {PhysicalSpecies::boron,
            Properties{PhysConst::m_p * amrex::Real(10.7319), PhysConst::q_e * amrex::Real(5.0)}},
        {PhysicalSpecies::boron10,
            Properties{PhysConst::m_p * amrex::Real(9.94060), PhysConst::q_e * amrex::Real(5.0)}},
        {PhysicalSpecies::boron11,
            Properties{PhysConst::m_p * amrex::Real(10.9298), PhysConst::q_e * amrex::Real(5.0)}},
        {PhysicalSpecies::carbon,
            Properties{PhysConst::m_e * amrex::Real(22032.0), PhysConst::q_e * amrex::Real(6.0)}},
        {PhysicalSpecies::nitrogen,
            Properties{PhysConst::m_e * amrex::Real(25716.9), PhysConst::q_e * amrex::Real(7.0)}},
        {PhysicalSpecies::oxygen,
            Properties{PhysConst::m_p * amrex::Real(15.8834), PhysConst::q_e * amrex::Real(8.0)}},
        {PhysicalSpecies::argon,
            Properties{PhysConst::m_p * amrex::Real(39.9480), PhysConst::q_e * amrex::Real(18.0)}},
        {PhysicalSpecies::copper,
            Properties{PhysConst::m_p * amrex::Real(63.0864), PhysConst::q_e * amrex::Real(29.0)}},
        {PhysicalSpecies::xenon,
            Properties{hysConst::m_p * amrex::Real(131.293), PhysConst::q_e * amrex::Real(54.0)}}
    };
}

namespace species
{
    std::optional<PhysicalSpecies> from_string(std::string species)
    {
        const auto phys_spec = string_to_species.find(species);
        return (phys_spec != species.end())?
            phys_spec : std::nullopt;
    }

    amrex::Real get_charge (PhysicalSpecies ps)
    {
        return species_to_properties.at(ps).charge;
    }

    amrex::Real get_mass (PhysicalSpecies ps)
    {
        return species_to_properties.at(ps).mass;
    }

    std::string get_name (PhysicalSpecies ps)
    {
        return species_to_string.at(ps);
    }
