/* Copyright 2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_SPECIESPHYSICALPROPERTIES_H_
#define WARPX_SPECIESPHYSICALPROPERTIES_H_

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_AmrCore.H>

#include <map>
#include <limits>

enum struct PhysicalSpecies{unspecified=0, electron, positron, photon, hydrogen, carbon, nitrogen, oxygen, copper};

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
}

#endif // WARPX_SPECIESPHYSICALPROPERTIES_H_
