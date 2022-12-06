/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef DSMC_FUNC_H_
#define DSMC_FUNC_H_

#include "Particles/Collision/ScatteringProcess.H"

#include <AMReX_DenseBins.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>


/**
 * \brief This functor perfoms DSMC (direct simulation Monte Carlo) collisions
 * within a cell. Particles are paired up and for each pair a stochastic process
 * determines whether a collision occurs. The algorith is the same as the one
 * used for binary Coulomb collisions and the nuclear fusion module.
 */
class DSMCFunc {
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the DSMCFunc class.
     */
    DSMCFunc () = default;

    /**
     * \brief Constructor of the DSMCFunc class
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc pointer to the MultiParticleContainer
     * @param[in] isSameSpecies whether the two colliding species are the same
     */
    DSMCFunc (const std::string collision_name,
              MultiParticleContainer const * const /*mypc*/,
              const bool isSameSpecies) : m_isSameSpecies(isSameSpecies)
    {
        using namespace amrex::literals;
        amrex::ParmParse pp_collision_name(collision_name);

        // query for a list of collision processes
        // these could be elastic, excitation, charge_exchange, back, etc.
        amrex::Vector<std::string> scattering_process_names;
        pp_collision_name.queryarr("scattering_processes", scattering_process_names);

        // create a vector of ScatteringProcess objects from each scattering
        // process name
        for (const auto& scattering_process : scattering_process_names) {
            std::string kw_cross_section = scattering_process + "_cross_section";
            std::string cross_section_file;
            pp_collision_name.query(kw_cross_section.c_str(), cross_section_file);

            amrex::ParticleReal energy = 0.0;
            // if the scattering process is excitation or ionization get the
            // energy associated with that process
            if (scattering_process.find("excitation") != std::string::npos ||
                scattering_process.find("ionization") != std::string::npos) {
                std::string kw_energy = scattering_process + "_energy";
                utils::parser::getWithParser(
                    pp_collision_name, kw_energy.c_str(), energy);
            }

            ScatteringProcess process(scattering_process, cross_section_file, energy);

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(process.type() != ScatteringProcessType::INVALID,
                                            "Cannot add an unknown scattering process type");

            m_scattering_processes.push_back(std::move(process));
        }
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

#if (defined WARPX_DIM_RZ)
            amrex::ParticleReal * const AMREX_RESTRICT theta1 = soa_1.m_rdata[PIdx::theta];
            amrex::ParticleReal * const AMREX_RESTRICT theta2 = soa_2.m_rdata[PIdx::theta];
#endif

            for (int k = 0; k < max_N; ++k)
            {
                // c1k : how many times the current particle of species 1 is paired with a particle
                // of species 2. Same for c2k.
                const int c1k = (k%NI1 < max_N%NI1) ? c1 + 1: c1;
                const int c2k = (k%NI2 < max_N%NI2) ? c2 + 1: c2;

#if (defined WARPX_DIM_RZ)
                /* In RZ geometry, macroparticles can collide with other macroparticles
                 * in the same *cylindrical* cell. For this reason, collisions between macroparticles
                 * are actually not local in space. In this case, the underlying assumption is that
                 * particles within the same cylindrical cell represent a cylindrically-symmetry
                 * momentum distribution function. Therefore, here, we temporarily rotate the
                 * momentum of one of the macroparticles in agreement with this cylindrical symmetry.
                 * (This is technically only valid if we use only the m=0 azimuthal mode in the simulation;
                 * there is a corresponding assert statement at initialization.) */
                amrex::ParticleReal const theta = theta2[I2[i2]]-theta1[I1[i1]];
                amrex::ParticleReal const u1xbuf = u1x[I1[i1]];
                u1x[I1[i1]] = u1xbuf*std::cos(theta) - u1y[I1[i1]]*std::sin(theta);
                u1y[I1[i1]] = u1xbuf*std::sin(theta) + u1y[I1[i1]]*std::cos(theta);
#endif

                // SingleNuclearFusionEvent(
                //     u1x[ I1[i1] ], u1y[ I1[i1] ], u1z[ I1[i1] ],
                //     u2x[ I2[i2] ], u2y[ I2[i2] ], u2z[ I2[i2] ],
                //     m1, m2, w1[ I1[i1] ]/c1k, w2[ I2[i2] ]/c2k,
                //     dt, dV, pair_index, p_mask, p_pair_reaction_weight,
                //     m_fusion_multiplier, multiplier_ratio,
                //     m_probability_threshold,
                //     m_probability_target_value,
                //     m_fusion_type, engine);

#if (defined WARPX_DIM_RZ)
                amrex::ParticleReal const u1xbuf_new = u1x[I1[i1]];
                u1x[I1[i1]] = u1xbuf_new*std::cos(-theta) - u1y[I1[i1]]*std::sin(-theta);
                u1y[I1[i1]] = u1xbuf_new*std::sin(-theta) + u1y[I1[i1]]*std::cos(-theta);
#endif

                p_pair_indices_1[pair_index] = I1[i1];
                p_pair_indices_2[pair_index] = I2[i2];
                ++i1; if ( i1 == static_cast<int>(I1e) ) { i1 = I1s; }
                ++i2; if ( i2 == static_cast<int>(I2e) ) { i2 = I2s; }
                ++pair_index;
            }

        }



private:
    amrex::Vector<ScatteringProcess> m_scattering_processes;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor> m_scattering_processes_exe;
    bool m_isSameSpecies;
};

#endif // DSMC_FUNC_H_