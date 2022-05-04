/* Copyright 2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_SPECIESPHYSICALPROPERTIES_H_
#define WARPX_SPECIESPHYSICALPROPERTIES_H_

#include <AMReX_REAL.H>

#include <optional>
#include <string>

enum struct PhysicalSpecies{unspecified=0, electron, positron, photon, hydrogen, helium, boron,
                            boron10, boron11, carbon, nitrogen, oxygen, argon, copper, xenon};

namespace species
{
    /**
     * \brief Returns the PhysicalSpecies associated to a given name
     *
     * \param[in] species_name the name of a species
     * \return the PhysicalSpecies corresponding to species_name (if it exists)
     */
    std::optional<PhysicalSpecies> from_string (const std::string& species_name);

    /**
     * \brief Returns the charge associated to a PhysicalSpecies
     *
     * \param[in] ps the PhysicalSpecies
     * \return the charge associated to the PhysicalSpecies
     */
    amrex::Real get_charge (const PhysicalSpecies& ps);

    /**
     * \brief Returns the mass associated to a PhysicalSpecies
     *
     * \param[in] ps the PhysicalSpecies
     * \return the mass associated to the PhysicalSpecies
     */
    amrex::Real get_mass (const PhysicalSpecies& ps);

    /**
     * \brief Returns the name associated to a PhysicalSpecies
     *
     * \param[in] ps the PhysicalSpecies
     * \return the name associated to the PhysicalSpecies
     */
    std::string get_name (const PhysicalSpecies& ps);
}

#endif // WARPX_SPECIESPHYSICALPROPERTIES_H_
