/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef NUCLEAR_FUSION_FUNC_H_
#define NUCLEAR_FUSION_FUNC_H_

#include "SingleNuclearFusionEvent.H"

#include "Particles/Collision/BinaryCollision/BinaryCollisionUtils.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_Algorithm.H>
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
     * @param[in] isSameSpecies whether the two colliding species are the same
     */
    NuclearFusionFunc (const std::string collision_name, MultiParticleContainer const * const mypc,
                       const bool isSameSpecies) : m_isSameSpecies(isSameSpecies)
    {
        using namespace amrex::literals;

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
        amrex::Abort("Nuclear fusion module does not currently work with single precision");
#endif

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
        // default fusion probability threshold
        m_probability_threshold = 0.02_rt;
        queryWithParser(pp_collision_name, "fusion_probability_threshold",
                                            m_probability_threshold);
        // default fusion probability target_value
        m_probability_target_value = 0.002_rt;
        queryWithParser(pp_collision_name, "fusion_probability_target_value",
                                            m_probability_target_value);
    }

    /**
     * \brief operator() of the NuclearFusionFunc class. Performs nuclear fusions at the cell level
     * using the algorithm described in Higginson et al., Journal of Computational Physics 388,
     * 439-453 (2019). Note that this function does not yet create the product particles, but
     * instead fills an array p_mask that stores which collisions result in a fusion event.
     *
     * Also note that there are three main differences between this implementation and the
     * algorithm described in Higginson's paper:
     * - 1) The transformation from the lab frame to the center of mass frame is nonrelativistic
     * in Higginson's paper. Here, we implement a relativistic generalization.
     * - 2) The behaviour when the estimated fusion probability is greater than one is not
     * specified in Higginson's paper. Here, we provide an implementation using two runtime
     * dependent parameters (fusion probability threshold and fusion probability target value). See
     * documentation for more details.
     * - 3) Here, we divide the weight of a particle by the number of times it is paired with other
     * particles. This was not explicitly specified in Higginson's paper.
     *
     * @param[in] I1s,I2s is the start index for I1,I2 (inclusive).
     * @param[in] I1e,I2e is the stop index for I1,I2 (exclusive).
     * @param[in] I1,I2 index arrays. They determine all elements that will be used.
     * @param[in] soa_1,soa_2 contain the struct of array data of the two species
     * @param[in] m1,m2 are masses.
     * @param[in] dt is the time step length between two collision calls.
     * @param[in] dV is the volume of the corresponding cell.
     * @param[in] cell_start_pair is the start index of the pairs in that cell.
     * @param[out] p_mask is a mask that will be set to true if a fusion event occurs for a given
     * pair. It is only needed here to store information that will be used later on when actually
     * creating the product particles.
     * @param[out] p_pair_indices_1,p_pair_indices_2 arrays that store the indices of the
     * particles of a given pair. They are only needed here to store information that will be used
     * later on when actually creating the product particles.
     * @param[out] p_pair_reaction_weight stores the weight of the product particles. It is only
     * needed here to store information that will be used later on when actually creating the
     * product particles.
     * @param[in] engine the random engine.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void operator() (
        index_type const I1s, index_type const I1e,
        index_type const I2s, index_type const I2e,
        index_type const* AMREX_RESTRICT I1,
        index_type const* AMREX_RESTRICT I2,
        SoaData_type soa_1, SoaData_type soa_2,
        GetParticlePosition /*get_position_1*/, GetParticlePosition /*get_position_2*/,
        amrex::ParticleReal const  /*q1*/, amrex::ParticleReal const  /*q2*/,
        amrex::ParticleReal const  m1, amrex::ParticleReal const  m2,
        amrex::Real const  dt, amrex::Real const dV,
        index_type const cell_start_pair, index_type* AMREX_RESTRICT p_mask,
        index_type* AMREX_RESTRICT p_pair_indices_1, index_type* AMREX_RESTRICT p_pair_indices_2,
        amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight,
        amrex::RandomEngine const& engine) const
        {

            amrex::ParticleReal * const AMREX_RESTRICT w1 = soa_1.m_rdata[PIdx::w];
            amrex::ParticleReal * const AMREX_RESTRICT u1x = soa_1.m_rdata[PIdx::ux];
            amrex::ParticleReal * const AMREX_RESTRICT u1y = soa_1.m_rdata[PIdx::uy];
            amrex::ParticleReal * const AMREX_RESTRICT u1z = soa_1.m_rdata[PIdx::uz];

            amrex::ParticleReal * const AMREX_RESTRICT w2 = soa_2.m_rdata[PIdx::w];
            amrex::ParticleReal * const AMREX_RESTRICT u2x = soa_2.m_rdata[PIdx::ux];
            amrex::ParticleReal * const AMREX_RESTRICT u2y = soa_2.m_rdata[PIdx::uy];
            amrex::ParticleReal * const AMREX_RESTRICT u2z = soa_2.m_rdata[PIdx::uz];

            // Number of macroparticles of each species
            const int NI1 = I1e - I1s;
            const int NI2 = I2e - I2s;
            const int max_N = amrex::max(NI1,NI2);

            int i1 = I1s;
            int i2 = I2s;
            int pair_index = cell_start_pair;

            // Because the number of particles of each species is not always equal (NI1 != NI2
            // in general), some macroparticles will be paired with multiple macroparticles of the
            // other species and we need to decrease their weight accordingly.
            // c1 corresponds to the minimum number of times a particle of species 1 will be paired
            // with a particle of species 2. Same for c2.
            const int c1 = amrex::max(NI2/NI1,1);
            const int c2 = amrex::max(NI1/NI2,1);

            // multiplier ratio to take into account unsampled pairs
            int multiplier_ratio;
            if (m_isSameSpecies)
            {
                multiplier_ratio = 2*max_N - 1;
            }
            else
            {
                multiplier_ratio = max_N;
            }

            for (int k = 0; k < max_N; ++k)
            {
                // c1k : how many times the current particle of species 1 is paired with a particle
                // of species 2. Same for c2k.
                const int c1k = (k%NI1 < max_N%NI1) ? c1 + 1: c1;
                const int c2k = (k%NI2 < max_N%NI2) ? c2 + 1: c2;
                SingleNuclearFusionEvent(
                    u1x[ I1[i1] ], u1y[ I1[i1] ], u1z[ I1[i1] ],
                    u2x[ I2[i2] ], u2y[ I2[i2] ], u2z[ I2[i2] ],
                    m1, m2, w1[ I1[i1] ]/c1k, w2[ I2[i2] ]/c2k,
                    dt, dV, pair_index, p_mask, p_pair_reaction_weight,
                    m_fusion_multiplier, multiplier_ratio,
                    m_probability_threshold,
                    m_probability_target_value,
                    m_fusion_type, engine);
                p_pair_indices_1[pair_index] = I1[i1];
                p_pair_indices_2[pair_index] = I2[i2];
                ++i1; if ( i1 == static_cast<int>(I1e) ) { i1 = I1s; }
                ++i2; if ( i2 == static_cast<int>(I2e) ) { i2 = I2s; }
                ++pair_index;
            }

        }

private:
    // Factor used to increase the number of fusion reaction by decreasing the weight of the
    // produced particles
    amrex::Real m_fusion_multiplier;
    // If the fusion multiplier is too high and results in a fusion probability that approaches
    // 1, there is a risk of underestimating the total fusion yield. In these cases, we reduce
    // the fusion multiplier used in a given collision. m_probability_threshold is the fusion
    // probability threshold above which we reduce the fusion multiplier.
    // m_probability_target_value is the target probability used to determine by how much
    // the fusion multiplier should be reduced.
    amrex::Real m_probability_threshold;
    amrex::Real m_probability_target_value;
    NuclearFusionType m_fusion_type;
    bool m_isSameSpecies;
};

#endif // NUCLEAR_FUSION_FUNC_H_
