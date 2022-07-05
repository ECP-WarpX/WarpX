/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "BinaryCollisionUtils.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <string>

namespace BinaryCollisionUtils{

    NuclearFusionType get_nuclear_fusion_type (const std::string collision_name,
                                               MultiParticleContainer const * const mypc)
        {
            amrex::ParmParse pp_collision_name(collision_name);
            amrex::Vector<std::string> species_names;
            pp_collision_name.getarr("species", species_names);
            auto& species1 = mypc->GetParticleContainerFromName(species_names[0]);
            auto& species2 = mypc->GetParticleContainerFromName(species_names[1]);

            if ((species1.AmIA<PhysicalSpecies::proton>() && species2.AmIA<PhysicalSpecies::boron11>())
                ||
                (species1.AmIA<PhysicalSpecies::boron11>() && species2.AmIA<PhysicalSpecies::proton>())
                )
            {
                return NuclearFusionType::ProtonBoron;
            }
            amrex::Abort("Binary nuclear fusion not implemented between species " +
                         species_names[0] + " of type " + species1.getSpeciesTypeName() +
                         " and species " + species_names[1] + " of type " +
                         species2.getSpeciesTypeName());
            return NuclearFusionType::Undefined;
        }

    CollisionType get_collision_type (const std::string collision_name,
                                      MultiParticleContainer const * const mypc)
        {
            amrex::ParmParse pp_collision_name(collision_name);
            std::string type;
            pp_collision_name.get("type", type);
            if (type == "nuclearfusion") {
                NuclearFusionType fusion_type = get_nuclear_fusion_type(collision_name, mypc);
                return nuclear_fusion_type_to_collision_type(fusion_type);
            }
            amrex::Abort(type + " is not a valid type of collision that creates new particles");
            return CollisionType::Undefined;
        }

    CollisionType nuclear_fusion_type_to_collision_type (const NuclearFusionType fusion_type)
        {
            if (fusion_type == NuclearFusionType::ProtonBoron)
                return CollisionType::ProtonBoronFusion;
            amrex::Abort("Invalid nuclear fusion type");
            return CollisionType::Undefined;
        }
}
