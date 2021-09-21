/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef NUCLEAR_FUSION_FUNC_H_
#define NUCLEAR_FUSION_FUNC_H_

#include "BinaryCollisionUtils.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_DenseBins.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * \brief This functor does binary nuclear fusions on a single cell.
 *  Particles of the two reacting species are paired with each other and for each pair we compute
 *  if a fusion event occurs. If so, we fill a mask (input parameter p_mask) with true so that
 *  product particles corresponding to a given pair can be effectively created in the particle
 *  creation functor.
 *  This functor also reads and contains the fusion multiplier.
 */
class NuclearFusionFunc{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the NuclearFusionFunc class.
     */
    NuclearFusionFunc () = default;

    /**
     * \brief Constructor of the NuclearFusionFunc class
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc pointer to the MultiParticleContainer
     */
    NuclearFusionFunc (const std::string collision_name, MultiParticleContainer const * const mypc)
    {
        using namespace amrex::literals;

        m_fusion_type = BinaryCollisionUtils::get_nuclear_fusion_type(collision_name, mypc);

        amrex::ParmParse pp_collision_name(collision_name);
        amrex::Vector<std::string> product_species_name;
        pp_collision_name.getarr("product_species", product_species_name);

        if (m_fusion_type == NuclearFusionType::ProtonBoron)
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                product_species_name.size() == 1,
                "ERROR: Proton-boron must contain exactly one product species");
            auto& product_species = mypc->GetParticleContainerFromName(product_species_name[0]);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                product_species.AmIA<PhysicalSpecies::helium>(),
                "ERROR: Product species of proton-boron fusion must be of type helium");
        }

        // default fusion multiplier
        m_fusion_multiplier = 1.0_rt;
        queryWithParser(pp_collision_name, "fusion_multiplier", m_fusion_multiplier);
    }

    /**
     * \brief operator() of the NuclearFusionFunc class. Needs to be implemented.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void operator() (
        index_type const /*I1s*/, index_type const /*I1e*/,
        index_type const /*I2s*/, index_type const /*I2e*/,
        index_type* /*I1*/,       index_type* /*I2*/,
        SoaData_type /*soa_1*/, SoaData_type /*soa_2*/,
        GetParticlePosition /*get_position_1*/, GetParticlePosition /*get_position_2*/,
        amrex::Real const  /*q1*/, amrex::Real const  /*q2*/,
        amrex::Real const  /*m1*/, amrex::Real const  /*m2*/,
        amrex::Real const  /*dt*/, amrex::Real const /*dV*/,
        index_type const /*cell_start_pair*/, index_type* /*p_mask*/,
        index_type* /*p_pair_indices_1*/, index_type* /*p_pair_indices_2*/,
        amrex::ParticleReal* /*p_pair_reaction_weight*/,
        amrex::RandomEngine const& /*engine*/) const
        {}

private:
    amrex::Real m_fusion_multiplier;
    NuclearFusionType m_fusion_type;
};

#endif // NUCLEAR_FUSION_FUNC_H_
