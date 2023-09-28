/* Copyright 2023 RemiLehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECIES_UTILS_H_
#define WARPX_SPECIES_UTILS_H_

#include <AMReX_REAL.H>
#include "Initialization/InjectorDensity.H"
#include "Initialization/InjectorMomentum.H"
#include "Particles/SpeciesPhysicalProperties.H"

namespace SpeciesUtils {

    void StringParseAbortMessage(const std::string& var,
                                 const std::string& name);

    void extractSpeciesProperties ( std::string const& species_name,
        std::string const& injection_style, amrex::Real& charge, amrex::Real& mass,
        PhysicalSpecies& physical_species);

    void parseDensity (std::string const& species_name,
        std::unique_ptr<InjectorDensity,InjectorDensityDeleter>& h_inj_rho);

    void parseMomentum (std::string const& species_name, const std::string& style,
        std::unique_ptr<InjectorMomentum,InjectorMomentumDeleter>& h_inj_mom,
        int flux_normal_axis=0, int flux_direction=0);

}

#endif // WARPX_SPECIES_UTILS_H_
