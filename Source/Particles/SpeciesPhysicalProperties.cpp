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
        {"muon"       , PhysicalSpecies::muon},
        {"antimuon"   , PhysicalSpecies::antimuon},
        {"photon"     , PhysicalSpecies::photon},
        {"hydrogen"   , PhysicalSpecies::hydrogen},
        {"protium"    , PhysicalSpecies::protium},
        {"deuterium"  , PhysicalSpecies::deuterium},
        {"tritium"    , PhysicalSpecies::tritium},
        {"proton"     , PhysicalSpecies::proton},
        {"helium"     , PhysicalSpecies::helium},
        {"alpha"      , PhysicalSpecies::helium},
        {"helium3"    , PhysicalSpecies::helium3},
        {"helium4"    , PhysicalSpecies::helium4},
        {"alpha"      , PhysicalSpecies::alpha},
        {"lithium"    , PhysicalSpecies::lithium},
        {"lithium6"   , PhysicalSpecies::lithium6},
        {"lithium7"   , PhysicalSpecies::lithium7},
        {"beryllium"  , PhysicalSpecies::beryllium},
        {"boron"      , PhysicalSpecies::boron},
        {"boron10"    , PhysicalSpecies::boron10},
        {"boron11"    , PhysicalSpecies::boron11},
        {"carbon"     , PhysicalSpecies::carbon},
        {"carbon12"   , PhysicalSpecies::carbon12},
        {"carbon13"   , PhysicalSpecies::carbon13},
        {"nitrogen"   , PhysicalSpecies::nitrogen},
        {"nitrogen14" , PhysicalSpecies::nitrogen14},
        {"nitrogen15" , PhysicalSpecies::nitrogen15},
        {"oxygen"     , PhysicalSpecies::oxygen},
        {"oxygen16"   , PhysicalSpecies::oxygen16},
        {"oxygen17"   , PhysicalSpecies::oxygen17},
        {"oxygen18"   , PhysicalSpecies::oxygen18},
        {"fluorine"   , PhysicalSpecies::fluorine},
        {"neon"       , PhysicalSpecies::neon},
        {"neon20"     , PhysicalSpecies::neon20},
        {"neon21"     , PhysicalSpecies::neon21},
        {"neon22"     , PhysicalSpecies::neon22},
        {"aluminium"  , PhysicalSpecies::aluminium},
        {"argon"      , PhysicalSpecies::argon},
        {"copper"     , PhysicalSpecies::copper},
        {"xenon"      , PhysicalSpecies::xenon},
        {"gold"       , PhysicalSpecies::gold},
    };

    const auto species_to_string = std::map<PhysicalSpecies, std::string>{
        {PhysicalSpecies::unspecified, "unspecified"},
        {PhysicalSpecies::electron    , "electron"},
        {PhysicalSpecies::positron   , "positron"},
        {PhysicalSpecies::muon       , "muon"},
        {PhysicalSpecies::antimuon   , "antimuon"},
        {PhysicalSpecies::photon     , "photon"},
        {PhysicalSpecies::hydrogen   , "hydrogen"},
        {PhysicalSpecies::protium    , "protium"},
        {PhysicalSpecies::deuterium  , "deuterium"},
        {PhysicalSpecies::tritium    , "tritium"},
        {PhysicalSpecies::proton     , "proton"},
        {PhysicalSpecies::helium     , "helium"},
        {PhysicalSpecies::helium3    , "helium3"},
        {PhysicalSpecies::helium4    , "helium4"},
        {PhysicalSpecies::alpha      , "alpha"},
        {PhysicalSpecies::lithium    , "lithium"},
        {PhysicalSpecies::lithium6   , "lithium6"},
        {PhysicalSpecies::lithium7   , "lithium7"},
        {PhysicalSpecies::beryllium  , "beryllium"},
        {PhysicalSpecies::boron      , "boron"},
        {PhysicalSpecies::boron10    , "boron10"},
        {PhysicalSpecies::boron11    , "boron11"},
        {PhysicalSpecies::carbon     , "carbon"},
        {PhysicalSpecies::carbon12   , "carbon12"},
        {PhysicalSpecies::carbon13   , "carbon13"},
        {PhysicalSpecies::nitrogen   , "nitrogen"},
        {PhysicalSpecies::nitrogen14 , "nitrogen14"},
        {PhysicalSpecies::nitrogen15 , "nitrogen15"},
        {PhysicalSpecies::oxygen     , "oxygen"},
        {PhysicalSpecies::fluorine   , "fluorine"},
        {PhysicalSpecies::neon       , "neon"},
        {PhysicalSpecies::neon20     , "neon20"},
        {PhysicalSpecies::neon21     , "neon21"},
        {PhysicalSpecies::neon22     , "neon22"},
        {PhysicalSpecies::aluminium  , "aluminium"},
        {PhysicalSpecies::argon      , "argon"},
        {PhysicalSpecies::copper     , "copper"},
        {PhysicalSpecies::xenon      , "xenon"},
        {PhysicalSpecies::gold       , "gold"}
    };

    constexpr auto quiet_NaN = std::numeric_limits<amrex::Real>::quiet_NaN();

    const
    std::map<PhysicalSpecies,Properties> species_to_properties
    {
        {PhysicalSpecies::unspecified,
            Properties{quiet_NaN, quiet_NaN}},
        {PhysicalSpecies::electron,
            Properties{                           PhysConst::m_e,                    -PhysConst::q_e}},
        {PhysicalSpecies::positron,
            Properties{                           PhysConst::m_e,                     PhysConst::q_e}},
        {PhysicalSpecies::muon,
            Properties{amrex::Real(206.7682830) * PhysConst::m_e,                    -PhysConst::q_e}},
        {PhysicalSpecies::antimuon,
            Properties{amrex::Real(206.7682830) * PhysConst::m_e,                     PhysConst::q_e}},
        {PhysicalSpecies::photon ,
            Properties{                         amrex::Real(0.0),                    amrex::Real(0.0)}},
        {PhysicalSpecies::hydrogen,
            Properties{      amrex::Real(1.008) * PhysConst::m_u,                     PhysConst::q_e}},
        {PhysicalSpecies::protium,
            Properties{     amrex::Real(1.0078) * PhysConst::m_u,                     PhysConst::q_e}},
        {PhysicalSpecies::deuterium,
            Properties{     amrex::Real(2.0141) * PhysConst::m_u,                     PhysConst::q_e}},
        {PhysicalSpecies::tritium,
            Properties{     amrex::Real(3.0160) * PhysConst::m_u,                     PhysConst::q_e}},
        {PhysicalSpecies::proton,
            Properties{                           PhysConst::m_p,                     PhysConst::q_e}},
        {PhysicalSpecies::helium,
            Properties{     amrex::Real(4.0026) * PhysConst::m_u,  amrex::Real(2.0) * PhysConst::q_e}},
        {PhysicalSpecies::helium3,
            Properties{     amrex::Real(3.0160) * PhysConst::m_u,  amrex::Real(2.0) * PhysConst::q_e}},
        {PhysicalSpecies::helium4,
            Properties{     amrex::Real(4.0026) * PhysConst::m_u,  amrex::Real(2.0) * PhysConst::q_e}},
        {PhysicalSpecies::alpha,
            Properties{     amrex::Real(4.0015) * PhysConst::m_u,  amrex::Real(2.0) * PhysConst::q_e}},
        {PhysicalSpecies::lithium,
            Properties{       amrex::Real(6.94) * PhysConst::m_u,  amrex::Real(3.0) * PhysConst::q_e}},
        {PhysicalSpecies::lithium6,
            Properties{     amrex::Real(6.0151) * PhysConst::m_u,  amrex::Real(3.0) * PhysConst::q_e}},
        {PhysicalSpecies::lithium7,
            Properties{     amrex::Real(7.0160) * PhysConst::m_u,  amrex::Real(3.0) * PhysConst::q_e}},
        {PhysicalSpecies::beryllium,
            Properties{     amrex::Real(9.0122) * PhysConst::m_u,  amrex::Real(4.0) * PhysConst::q_e}},
        {PhysicalSpecies::boron,
            Properties{      amrex::Real(10.81) * PhysConst::m_u,  amrex::Real(5.0) * PhysConst::q_e}},
        {PhysicalSpecies::boron10,
            Properties{    amrex::Real(10.0129) * PhysConst::m_u,  amrex::Real(5.0) * PhysConst::q_e}},
        {PhysicalSpecies::boron11,
            Properties{    amrex::Real(11.0093) * PhysConst::m_u,  amrex::Real(5.0) * PhysConst::q_e}},
        {PhysicalSpecies::carbon,
            Properties{     amrex::Real(12.011) * PhysConst::m_u,  amrex::Real(6.0) * PhysConst::q_e}},
        {PhysicalSpecies::carbon12,
            Properties{         amrex::Real(12) * PhysConst::m_u,  amrex::Real(6.0) * PhysConst::q_e}},
        {PhysicalSpecies::carbon13,
            Properties{    amrex::Real(13.0033) * PhysConst::m_u,  amrex::Real(6.0) * PhysConst::q_e}},
        {PhysicalSpecies::nitrogen,
            Properties{     amrex::Real(14.007) * PhysConst::m_u,  amrex::Real(7.0) * PhysConst::q_e}},
        {PhysicalSpecies::nitrogen14,
            Properties{    amrex::Real(14.0031) * PhysConst::m_u,  amrex::Real(7.0) * PhysConst::q_e}},
        {PhysicalSpecies::nitrogen15,
            Properties{    amrex::Real(15.0001) * PhysConst::m_u,  amrex::Real(7.0) * PhysConst::q_e}},
        {PhysicalSpecies::oxygen,
            Properties{     amrex::Real(15.999) * PhysConst::m_u,  amrex::Real(8.0) * PhysConst::q_e}},
        {PhysicalSpecies::oxygen16,
            Properties{    amrex::Real(15.9949) * PhysConst::m_u,  amrex::Real(8.0) * PhysConst::q_e}},
        {PhysicalSpecies::oxygen17,
            Properties{    amrex::Real(16.9991) * PhysConst::m_u,  amrex::Real(8.0) * PhysConst::q_e}},
        {PhysicalSpecies::oxygen18,
            Properties{    amrex::Real(17.9992) * PhysConst::m_u,  amrex::Real(8.0) * PhysConst::q_e}},
        {PhysicalSpecies::fluorine,
            Properties{    amrex::Real(18.9984) * PhysConst::m_u,  amrex::Real(9.0) * PhysConst::q_e}},
        {PhysicalSpecies::neon,
            Properties{    amrex::Real(20.1798) * PhysConst::m_u, amrex::Real(10.0) * PhysConst::q_e}},
        {PhysicalSpecies::neon20,
            Properties{    amrex::Real(19.9924) * PhysConst::m_u, amrex::Real(10.0) * PhysConst::q_e}},
        {PhysicalSpecies::neon21,
            Properties{    amrex::Real(20.9938) * PhysConst::m_u, amrex::Real(10.0) * PhysConst::q_e}},
        {PhysicalSpecies::neon22,
            Properties{    amrex::Real(21.9914) * PhysConst::m_u, amrex::Real(10.0) * PhysConst::q_e}},
        {PhysicalSpecies::aluminium,
            Properties{    amrex::Real(26.9815) * PhysConst::m_u, amrex::Real(13.0) * PhysConst::q_e}},
        {PhysicalSpecies::argon,
            Properties{    amrex::Real(39.7924) * PhysConst::m_u, amrex::Real(18.0) * PhysConst::q_e}},
        {PhysicalSpecies::copper,
            Properties{     amrex::Real(63.546) * PhysConst::m_u, amrex::Real(29.0) * PhysConst::q_e}},
        {PhysicalSpecies::xenon,
            Properties{    amrex::Real(131.293) * PhysConst::m_u, amrex::Real(54.0) * PhysConst::q_e}},
        {PhysicalSpecies::gold,
            Properties{   amrex::Real(196.9666) * PhysConst::m_u, amrex::Real(79.0) * PhysConst::q_e}}
    };
}

namespace species
{
    std::optional<PhysicalSpecies> from_string(const std::string& species)
    {
        const auto phys_spec = string_to_species.find(species);
        return (phys_spec != string_to_species.end())?
            std::make_optional(phys_spec->second) : std::nullopt;
    }

    amrex::Real get_charge (const PhysicalSpecies& ps)
    {
        return species_to_properties.at(ps).charge;
    }

    amrex::Real get_mass (const PhysicalSpecies& ps)
    {
        return species_to_properties.at(ps).mass;
    }

    std::string get_name (const PhysicalSpecies& ps)
    {
        return species_to_string.at(ps);
    }
}
