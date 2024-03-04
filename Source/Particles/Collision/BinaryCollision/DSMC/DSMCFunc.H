/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef DSMC_FUNC_H_
#define DSMC_FUNC_H_

#include "CollisionFilterFunc.H"

#include "Particles/Collision/BinaryCollision/BinaryCollisionUtils.H"
#include "Particles/Collision/BinaryCollision/ShuffleFisherYates.H"
#include "Particles/Collision/CollisionBase.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX_DenseBins.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>

/**
 * \brief This class performs DSMC (direct simulation Monte Carlo) collisions
 * within a cell. Particles are paired up and for each pair a stochastic process
 * determines whether a collision occurs. The algorithm is similar to the one
 * used for binary Coulomb collisions and the nuclear fusion module.
 */
class DSMCFunc
{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = ParticleTileType::ParticleTileDataType;
    using ParticleBins = amrex::DenseBins<ParticleTileDataType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:

    /**
     * \brief Constructor of the DSMCFunc class.
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
              [[maybe_unused]] MultiParticleContainer const * const mypc,
              [[maybe_unused]] const bool isSameSpecies)
        //   : m_isSameSpecies{isSameSpecies}
    {
        using namespace amrex::literals;

        const amrex::ParmParse pp_collision_name(collision_name);

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

            // if the scattering process is excitation or ionization get the
            // energy associated with that process
            amrex::ParticleReal energy = 0._prt;
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

#ifdef AMREX_USE_GPU
        amrex::Gpu::HostVector<ScatteringProcess::Executor> h_scattering_processes_exe;
        for (auto const& p : m_scattering_processes) {
            h_scattering_processes_exe.push_back(p.executor());
        }
        m_scattering_processes_exe.resize(h_scattering_processes_exe.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_scattering_processes_exe.begin(),
                            h_scattering_processes_exe.end(), m_scattering_processes_exe.begin());
        amrex::Gpu::streamSynchronize();
#else
        for (auto const& p : m_scattering_processes) {
            m_scattering_processes_exe.push_back(p.executor());
        }
#endif
        m_process_count = static_cast<int>(m_scattering_processes_exe.size());
    }

    /**
     * \brief operator() of the DSMCFunc class. Performs DSMC collisions at the cell level.
     * Note that this function does not yet create the product particles, but
     * instead fills an array p_mask that stores which pairs result in a collision event.
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
        GetParticlePosition<PIdx> /*get_position_1*/, GetParticlePosition<PIdx> /*get_position_2*/,
        amrex::ParticleReal const  /*q1*/, amrex::ParticleReal const  /*q2*/,
        amrex::ParticleReal const  m1, amrex::ParticleReal const  m2,
        amrex::Real const  dt, amrex::Real const dV,
        index_type const cell_start_pair, index_type* AMREX_RESTRICT p_mask,
        index_type* AMREX_RESTRICT p_pair_indices_1, index_type* AMREX_RESTRICT p_pair_indices_2,
        amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight,
        amrex::RandomEngine const& engine) const
    {
        auto* scattering_processes = m_scattering_processes_exe.data();

        amrex::ParticleReal * const AMREX_RESTRICT w1 = soa_1.m_rdata[PIdx::w];
        amrex::ParticleReal * const AMREX_RESTRICT u1x = soa_1.m_rdata[PIdx::ux];
        amrex::ParticleReal * const AMREX_RESTRICT u1y = soa_1.m_rdata[PIdx::uy];
        amrex::ParticleReal * const AMREX_RESTRICT u1z = soa_1.m_rdata[PIdx::uz];

        amrex::ParticleReal * const AMREX_RESTRICT w2 = soa_2.m_rdata[PIdx::w];
        amrex::ParticleReal * const AMREX_RESTRICT u2x = soa_2.m_rdata[PIdx::ux];
        amrex::ParticleReal * const AMREX_RESTRICT u2y = soa_2.m_rdata[PIdx::uy];
        amrex::ParticleReal * const AMREX_RESTRICT u2z = soa_2.m_rdata[PIdx::uz];

        // Number of macroparticles of each species
        const index_type NI1 = I1e - I1s;
        const index_type NI2 = I2e - I2s;
        const index_type max_N = amrex::max(NI1,NI2);

        index_type i1 = I1s;
        index_type i2 = I2s;
        index_type pair_index = cell_start_pair;

        // Because the number of particles of each species is not always equal (NI1 != NI2
        // in general), some macroparticles will be paired with multiple macroparticles of the
        // other species and we need to decrease their weight accordingly.
        // c1 corresponds to the minimum number of times a particle of species 1 will be paired
        // with a particle of species 2. Same for c2.
        // index_type(1): https://github.com/AMReX-Codes/amrex/pull/3684
        const index_type c1 = amrex::max(NI2/NI1, index_type(1));
        const index_type c2 = amrex::max(NI1/NI2, index_type(1));

#if (defined WARPX_DIM_RZ)
        amrex::ParticleReal * const AMREX_RESTRICT theta1 = soa_1.m_rdata[PIdx::theta];
        amrex::ParticleReal * const AMREX_RESTRICT theta2 = soa_2.m_rdata[PIdx::theta];
#endif

        for (index_type k = 0; k < max_N; ++k)
        {
            // c1k : how many times the current particle of species 1 is paired with a particle
            // of species 2. Same for c2k.
            const index_type c1k = (k%NI1 < max_N%NI1) ? c1 + 1: c1;
            const index_type c2k = (k%NI2 < max_N%NI2) ? c2 + 1: c2;

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

            CollisionPairFilter(
                u1x[ I1[i1] ], u1y[ I1[i1] ], u1z[ I1[i1] ],
                u2x[ I2[i2] ], u2y[ I2[i2] ], u2z[ I2[i2] ],
                m1, m2, w1[ I1[i1] ]/c1k, w2[ I2[i2] ]/c2k,
                dt, dV, pair_index, p_mask, p_pair_reaction_weight,
                max_N, m_process_count, scattering_processes, engine);

#if (defined WARPX_DIM_RZ)
            amrex::ParticleReal const u1xbuf_new = u1x[I1[i1]];
            u1x[I1[i1]] = u1xbuf_new*std::cos(-theta) - u1y[I1[i1]]*std::sin(-theta);
            u1y[I1[i1]] = u1xbuf_new*std::sin(-theta) + u1y[I1[i1]]*std::cos(-theta);
#endif

            p_pair_indices_1[pair_index] = I1[i1];
            p_pair_indices_2[pair_index] = I2[i2];
            ++i1; if ( i1 == I1e ) { i1 = I1s; }
            ++i2; if ( i2 == I2e ) { i2 = I2s; }
            ++pair_index;
        }
    }

private:
    amrex::Vector<ScatteringProcess> m_scattering_processes;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor> m_scattering_processes_exe;
    int m_process_count;
    // bool m_isSameSpecies;
};

#endif // DSMC_FUNC_H_
