/* Copyright 2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_SPECIESPHYSICALPROPERTIES_H_
#define WARPX_SPECIESPHYSICALPROPERTIES_H_

#include "Utils/WarpXConst.H"

#include <AMReX_AmrCore.H>
#include <AMReX_REAL.H>

#include <limits>
#include <map>
#include <string>

enum struct PhysicalSpecies{unspecified=0, electron, positron, photon, hydrogen, helium, boron,
                            carbon, nitrogen, oxygen, copper};

namespace species
{
    AMREX_FORCE_INLINE
    PhysicalSpecies from_string(std::string species)
    {
        if( species=="unspecified" )
            return PhysicalSpecies::unspecified;
        if( species=="electron" )
            return PhysicalSpecies::electron;
        if( species=="positron" )
            return PhysicalSpecies::positron;
        if( species=="photon" )
            return PhysicalSpecies::photon;
        if( species=="hydrogen" )
            return PhysicalSpecies::hydrogen;
        if( species=="proton" )
            return PhysicalSpecies::hydrogen;
        if( species=="helium" )
            return PhysicalSpecies::helium;
        if( species=="alpha" )
            return PhysicalSpecies::helium;
        if( species=="boron" )
            return PhysicalSpecies::boron;
        if( species=="carbon" )
            return PhysicalSpecies::carbon;
        if( species=="nitrogen" )
            return PhysicalSpecies::nitrogen;
        if( species=="oxygen" )
            return PhysicalSpecies::oxygen;
        if( species=="copper" )
            return PhysicalSpecies::copper;
        amrex::Abort("unknown PhysicalSpecies");
        return PhysicalSpecies::unspecified;
    }

    AMREX_FORCE_INLINE
    amrex::Real get_charge (PhysicalSpecies ps)
    {
        switch(ps) {
        case PhysicalSpecies::unspecified:
            return std::numeric_limits<amrex::Real>::quiet_NaN();
        case PhysicalSpecies::electron:
            return -PhysConst::q_e;
        case PhysicalSpecies::positron:
            return PhysConst::q_e;
        case PhysicalSpecies::photon:
            return 0.;
        case PhysicalSpecies::hydrogen:
            return PhysConst::q_e;
        case PhysicalSpecies::helium:
            return PhysConst::q_e * amrex::Real(2.0);
        case PhysicalSpecies::boron:
            return PhysConst::q_e * amrex::Real(5.0);
        case PhysicalSpecies::carbon:
            return PhysConst::q_e * amrex::Real(6.0);
        case PhysicalSpecies::nitrogen:
            return PhysConst::q_e * amrex::Real(7.0);
        case PhysicalSpecies::oxygen:
            return PhysConst::q_e * amrex::Real(8.0);
        case PhysicalSpecies::copper:
            return PhysConst::q_e * amrex::Real(29.0);
        default:
            amrex::Abort("unknown PhysicalSpecies");
            return 0.;
        }
    }

    AMREX_FORCE_INLINE
    amrex::Real get_mass (PhysicalSpecies ps)
    {
        switch(ps) {
        case PhysicalSpecies::unspecified:
            return std::numeric_limits<amrex::Real>::quiet_NaN();
        case PhysicalSpecies::electron:
            return PhysConst::m_e;
        case PhysicalSpecies::positron:
            return PhysConst::m_e;
        case PhysicalSpecies::photon:
            return 0.;
        case PhysicalSpecies::hydrogen:
            return PhysConst::m_p;
        case PhysicalSpecies::helium:
            return PhysConst::m_p * amrex::Real(3.97369);
        case PhysicalSpecies::boron:
            return PhysConst::m_p * amrex::Real(10.7319);
        case PhysicalSpecies::carbon:
            return PhysConst::m_e * amrex::Real(22032.0);
        case PhysicalSpecies::nitrogen:
            return PhysConst::m_e * amrex::Real(25716.9);
        case PhysicalSpecies::oxygen:
            return PhysConst::m_p * amrex::Real(15.8834);
        case PhysicalSpecies::copper:
            return PhysConst::m_p * amrex::Real(63.0864);
        default:
            amrex::Abort("unknown PhysicalSpecies");
            return 0.;
        }
    }

    AMREX_FORCE_INLINE
    std::string get_name (PhysicalSpecies ps)
    {
        switch(ps) {
        case PhysicalSpecies::unspecified:
            return "unspecified";
        case PhysicalSpecies::electron:
            return "electron";
        case PhysicalSpecies::positron:
            return "positron";
        case PhysicalSpecies::photon:
            return "photon";
        case PhysicalSpecies::hydrogen:
            return "hydrogen";
        case PhysicalSpecies::helium:
            return "helium";
        case PhysicalSpecies::boron:
            return "boron";
        case PhysicalSpecies::carbon:
            return "carbon";
        case PhysicalSpecies::nitrogen:
            return "nitrogen";
        case PhysicalSpecies::oxygen:
            return "oxygen";
        case PhysicalSpecies::copper:
            return "copper";
        default:
            amrex::Abort("unknown PhysicalSpecies");
            return "";
        }
    }
}

#endif // WARPX_SPECIESPHYSICALPROPERTIES_H_
