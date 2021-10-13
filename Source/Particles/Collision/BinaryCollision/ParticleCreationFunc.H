/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PARTICLE_CREATION_FUNC_H_
#define PARTICLE_CREATION_FUNC_H_

#include "BinaryCollisionUtils.H"

#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_DenseBins.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * \brief This functor creates particles produced from a binary collision and sets their initial
 * properties (position, momentum, weight).
 */
class ParticleCreationFunc{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the ParticleCreationFunc class.
     */
    ParticleCreationFunc () = default;

    /**
     * \brief Constructor of the ParticleCreationFunc class
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc pointer to the MultiParticleContainer
     */
    ParticleCreationFunc (const std::string collision_name, MultiParticleContainer const * const mypc);

    /**
     * \brief operator() of the ParticleCreationFunc class. Needs to be implemented.
     */
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
                    ) const
        {
            return amrex::Vector<int>(m_num_product_species, 0);
        }

private:
    // How many different type of species the collision produces
    int m_num_product_species;
    // Vector of size m_num_product_species storing how many particles of a given species are
    // produced by a collision event
    amrex::Gpu::DeviceVector<int> m_num_products;
    CollisionType m_collision_type;
};


/**
 * \brief This class does nothing and is used as second template parameter for binary collisions
 * that do not create particles.
 */
class NoParticleCreationFunc{
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    NoParticleCreationFunc () = default;

    NoParticleCreationFunc (const std::string /*collision_name*/,
                            MultiParticleContainer const * const /*mypc*/) {}

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
                    ) const
        {
            return amrex::Vector<int>();
        }
};

#endif // PARTICLE_CREATION_FUNC_H_
