/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef BINARY_COLLISION_UTILS_H_
#define BINARY_COLLISION_UTILS_H_

#include <string>

#include "Particles/MultiParticleContainer.H"

enum struct CollisionType { DeuteriumTritiumToNeutronHeliumFusion,
                            DeuteriumDeuteriumToProtonTritiumFusion,
                            DeuteriumDeuteriumToNeutronHeliumFusion,
                            DeuteriumHeliumToProtonHeliumFusion,
                            ProtonBoronToAlphasFusion,
                            Undefined };

enum struct NuclearFusionType {
                                DeuteriumTritiumToNeutronHelium,
                                DeuteriumDeuteriumToProtonTritium,
                                DeuteriumDeuteriumToNeutronHelium,
                                DeuteriumHeliumToProtonHelium,
                                ProtonBoronToAlphas,
                                Undefined };

namespace BinaryCollisionUtils{

    NuclearFusionType get_nuclear_fusion_type (std::string collision_name,
                                               MultiParticleContainer const * mypc);

    CollisionType get_collision_type (std::string collision_name,
                                      MultiParticleContainer const * mypc);

    CollisionType nuclear_fusion_type_to_collision_type (NuclearFusionType fusion_type);

    void get_collision_parameters (
        const amrex::ParticleReal& u1x, const amrex::ParticleReal& u1y,
        const amrex::ParticleReal& u1z, const amrex::ParticleReal& u2x,
        const amrex::ParticleReal& u2y, const amrex::ParticleReal& u2z,
        const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
        amrex::ParticleReal& E_coll, amrex::ParticleReal& v_coll,
        amrex::ParticleReal& lab_to_COM_factor
    );
}

#endif // BINARY_COLLISION_UTILS_H_
