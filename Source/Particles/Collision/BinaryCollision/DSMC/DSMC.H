/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef DSMC_H_
#define DSMC_H_

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
 * \brief This class perfoms DSMC (direct simulation Monte Carlo) collisions
 * within a cell. Particles are paired up and for each pair a stochastic process
 * determines whether a collision occurs. The algorithm is similar to the one
 * used for binary Coulomb collisions and the nuclear fusion module.
 */
class DSMC final
    : public CollisionBase
{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;
    using index_type = ParticleBins::index_type;

public:
    /**
     * \brief Default constructor of the DSMCFunc class.
     */
    virtual ~DSMC () = default;

    /**
     * \brief Constructor of the DSMCFunc class
     *
     * @param[in] collision_name the name of the collision
     */
    DSMC (const std::string collision_name);

    /** Perform the collisions
     *
     * @param cur_time Current time
     * @param dt Time step size
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc) override;

    /** Perform all binary collisions within a tile
     *
     * \param[in] lev the mesh-refinement level
     * \param[in] mfi iterator for multifab
     * \param species_1 first species container
     * \param species_2 second species container
     * \param copy_species1 SmartCopy for species_1
     * \param copy_species2 SmartCopy for species_2
     *
     */
    void doCollisionsWithinTile (
        amrex::Real dt, int const lev, amrex::MFIter const& mfi,
        WarpXParticleContainer& species_1,
        WarpXParticleContainer& species_2,
        SmartCopy& copy_species1,
        SmartCopy& copy_species2 );

    /**
     * \brief This function determines whether a collision occurs for a given
     * pair of particles.
     *
     * @param[in] u1x,u1y,u1z momenta of the first colliding particle
     * @param[in] u2x,u2y,u2z momenta of the second colliding particle
     * @param[in] m1,m2 masses
     * @param[in] w1,w2 effective weight of the colliding particles
     * @param[in] dt is the time step length between two collision calls.
     * @param[in] dV is the volume of the corresponding cell.
     * @param[in] pair_index is the index of the colliding pair
     * @param[out] p_mask is a mask that will be set to true if fusion occurs for that pair
     * @param[out] p_pair_reaction_weight stores the weight of the product particles
     * @param[in] engine the random engine.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void collision_pair_filter (const amrex::ParticleReal& u1x, const amrex::ParticleReal& u1y,
                                const amrex::ParticleReal& u1z, const amrex::ParticleReal& u2x,
                                const amrex::ParticleReal& u2y, const amrex::ParticleReal& u2z,
                                const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
                                amrex::ParticleReal w1, amrex::ParticleReal w2,
                                const amrex::Real& dt, const amrex::ParticleReal& dV, const int& pair_index,
                                index_type* AMREX_RESTRICT p_mask,
                                amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight,
                                const int& multiplier_ratio,
                                const amrex::RandomEngine& engine) const
    {
        amrex::ParticleReal E_coll, v_coll, lab_to_COM_factor;

        const amrex::ParticleReal w_min = amrex::min(w1, w2);
        const amrex::ParticleReal w_max = amrex::max(w1, w2);

        BinaryCollisionUtils::get_collision_parameters(
            u1x, u1y, u1z, u2x, u2y, u2z, m1, m2,
            E_coll, v_coll, lab_to_COM_factor);

        // convert E_coll to eV
        E_coll /= PhysConst::q_e;

        amrex::ParticleReal sigma_tot = 0._prt;
        for (int ii = 0; ii < m_scattering_processes_exe.size(); ii++) {
            sigma_tot += m_scattering_processes_exe.data()[ii].getCrossSection(E_coll);
        }

        // calculate total collision probability
        amrex::ParticleReal exponent = (
            lab_to_COM_factor * multiplier_ratio * w_max
            * sigma_tot * v_coll * dt / dV
        );

        // Compute actual collision probability that is always between zero and one
        // In principle this is obtained by computing 1 - exp(-probability_estimate)
        // However, the computation of this quantity can fail numerically when probability_estimate is
        // too small (e.g. exp(-probability_estimate) returns 1 and the computation returns 0).
        // In this case, we simply use "probability_estimate" instead of 1 - exp(-probability_estimate)
        // The threshold exp_threshold at which we switch between the two formulas is determined by the
        // fact that computing the exponential is only useful if it can resolve the x^2/2 term of its
        // Taylor expansion, i.e. the square of probability_estimate should be greater than the
        // machine epsilon.
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
        constexpr auto exp_threshold = amrex::ParticleReal(1.e-3);
#else
        constexpr auto exp_threshold = amrex::ParticleReal(5.e-8);
#endif
        const amrex::ParticleReal probability = (exponent < exp_threshold) ?
                                    exponent: 1._prt - std::exp(-exponent);

        // Now we determine if a collision should occur
        if (amrex::Random(engine) < probability)
        {
            const amrex::ParticleReal random_number = amrex::Random(engine);
            amrex::ParticleReal sigma = 0._prt;
            for (int ii = 0; ii < m_scattering_processes_exe.size(); ii++) {
                sigma += m_scattering_processes_exe.data()[ii].getCrossSection(E_coll);
                if (random_number <= sigma / sigma_tot)
                {
                    p_mask[pair_index] = int(m_scattering_processes[ii].type());
                    p_pair_reaction_weight[pair_index] = w_min;
                    break;
                }
            }
        }
        else
        {
            p_mask[pair_index] = false;
        }
    }

    /**
     * \brief Function that determines if a collision occurs and if so, what
     * type.
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
    void collision_filter (
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

                collision_pair_filter(
                    u1x[ I1[i1] ], u1y[ I1[i1] ], u1z[ I1[i1] ],
                    u2x[ I2[i2] ], u2y[ I2[i2] ], u2z[ I2[i2] ],
                    m1, m2, w1[ I1[i1] ]/c1k, w2[ I2[i2] ]/c2k,
                    dt, dV, pair_index, p_mask, p_pair_reaction_weight,
                    max_N, engine);

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
};

#endif // DSMC_H_