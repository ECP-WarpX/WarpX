/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef BINARY_COLLISION_UTILS_H_
#define BINARY_COLLISION_UTILS_H_

#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_DenseBins.H>
#include <AMReX_Vector.H>

/**
 * \brief This class does nothing and is used as second template parameter for binary collisions
 * that do not create particles.
 */
class NoParticleCreationFunc{
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;

public:
    AMREX_INLINE
    amrex::Vector<int> operator() (
                    const SoaData_type /*soa_1*/, const SoaData_type /*soa_2*/,
                    const SoaData_type* /*soa_products*/,
                    const GetParticlePosition /*get_position_1*/,
                    const GetParticlePosition /*get_position_2*/,
                    const GetParticlePosition* /*get_position_products*/,
                    const index_type* /*p_mask*/, const index_type* /*products_np*/,
                    const SmartCopy* /*copy_species1*/, const SmartCopy* /*copy_species2*/,
                    const index_type* /*p_pair_indices_1*/, const index_type* /*p_pair_indices_2*/,
                    const amrex::ParticleReal* /*p_pair_reaction_weight*/
                    ) const {
                        return amrex::Vector<int>();
                    }
};

#endif // BINARY_COLLISION_UTILS_H_
