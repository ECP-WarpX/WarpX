/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PROTON_BORON_FUSION_FUNC_H_
#define PROTON_BORON_FUSION_FUNC_H_

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_DenseBins.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>


/**
 * \brief This functor performs proton boron fusion on a single cell.
 *  It also reads and contains the fusion multiplier.
 */
class ProtonBoronFusionFunc{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the ProtonBoronFusionFunc class.
     */
    ProtonBoronFusionFunc () = default;

    /**
     * \brief Constructor of the ProtonBoronFusionFunc class
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc pointer to the MultiParticleContainer
     */
    ProtonBoronFusionFunc (const std::string collision_name, MultiParticleContainer* const mypc)
    {
        using namespace amrex::literals;
        amrex::ParmParse pp_collision_name(collision_name);

        amrex::Vector<std::string> species_names;
        pp_collision_name.getarr("species", species_names);
        auto& species1 = mypc->GetParticleContainerFromName(species_names[0]);
        auto& species2 = mypc->GetParticleContainerFromName(species_names[1]);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                species1.AmIA<PhysicalSpecies::hydrogen>(),
                "ERROR: First species of proton-boron fusion must be of type hydrogen");

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                species2.AmIA<PhysicalSpecies::boron>(),
                "ERROR: Second species of proton-boron fusion must be of type boron");

        std::string product_species_name;
        pp_collision_name.get("product_species", product_species_name);
        auto& product_species = mypc->GetParticleContainerFromName(product_species_name);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                product_species.AmIA<PhysicalSpecies::helium>(),
                "ERROR: Product species of proton-boron fusion must be of type helium");

        // default fusion multiplier
        m_fusion_multiplier = 1.0_rt;
        queryWithParser(pp_collision_name, "fusion_multiplier", m_fusion_multiplier);
    }

    /**
     * \brief operator() of the ProtonBoronFusionFunc class. Needs to be implemented.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void operator() (
        index_type const /*I1s*/, index_type const /*I1e*/,
        index_type const /*I2s*/, index_type const /*I2e*/,
        index_type* /*I1*/,       index_type* /*I2*/,
        SoaData_type /*soa_1*/, SoaData_type /*soa_2*/, SoaData_type /*soa_product*/,
        GetParticlePosition /*get_position_1*/, GetParticlePosition /*get_position_2*/,
        GetParticlePosition /*get_position_product*/,
        amrex::Real const  /*q1*/, amrex::Real const  /*q2*/,
        amrex::Real const  /*m1*/, amrex::Real const  /*m2*/,
        amrex::Real const  /*dt*/, amrex::Real const /*dV*/,
        amrex::RandomEngine const& /*engine*/) const
        {}

private:
    amrex::Real m_fusion_multiplier;
};

#endif // PROTON_BORON_FUSION_FUNC_H_
