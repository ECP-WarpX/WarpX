/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PAIRWISE_COULOMB_COLLISION_FUNC_H_
#define PAIRWISE_COULOMB_COLLISION_FUNC_H_

#include "ElasticCollisionPerez.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXUtil.H"

#include <AMReX_DenseBins.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>


/**
 * \brief This functor performs pairwise Coulomb collision on a single cell by calling the function
 *        ElasticCollisionPerez. It also reads and contains the Coulomb logarithm.
 */
class PairWiseCoulombCollisionFunc{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the PairWiseCoulombCollisionFunc class.
     */
    PairWiseCoulombCollisionFunc () = default;

    /**
     * \brief Constructor of the PairWiseCoulombCollisionFunc class
     *
     * @param[in] collision_name the name of the collision
     */
    PairWiseCoulombCollisionFunc (const std::string collision_name,
                                  MultiParticleContainer const * const /*mypc*/,
                                  const bool /*isSameSpecies*/)
    {
        using namespace amrex::literals;
        amrex::ParmParse pp_collision_name(collision_name);
        // default Coulomb log, if < 0, will be computed automatically
        m_CoulombLog = -1.0_rt;
        queryWithParser(pp_collision_name, "CoulombLog", m_CoulombLog);
    }

    /**
     * \brief operator() of the PairWiseCoulombCollisionFunc class. Performs Coulomb collisions
     * at the cell level by calling ElasticCollisionPerez.
     *
     * @param[in] I1s,I2s is the start index for I1,I2 (inclusive).
     * @param[in] I1e,I2e is the stop index for I1,I2 (exclusive).
     * @param[in] I1 and I2 are the index arrays. They determine all elements that will be used.
     * @param[in, out] soa_1 and soa_2 contain the struct of array data of the two species.
     * @param[in] q1 and q2 are charges. m1 and m2 are masses.
     * @param[in] dt is the time step length between two collision calls.
     * @param[in] dV is the volume of the corresponding cell.
     * @param[in] engine the random engine.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void operator() (
        index_type const I1s, index_type const I1e,
        index_type const I2s, index_type const I2e,
        index_type* I1,       index_type* I2,
        SoaData_type soa_1, SoaData_type soa_2,
        GetParticlePosition /*get_position_1*/, GetParticlePosition /*get_position_2*/,
        amrex::Real const  q1, amrex::Real const  q2,
        amrex::Real const  m1, amrex::Real const  m2,
        amrex::Real const  dt, amrex::Real const dV,
        index_type const /*cell_start_pair*/, index_type* /*p_mask*/,
        index_type* /*p_pair_indices_1*/, index_type* /*p_pair_indices_2*/,
        amrex::ParticleReal* /*p_pair_reaction_weight*/,
        amrex::RandomEngine const& engine) const
        {
            ElasticCollisionPerez(
                    I1s, I1e, I2s, I2e, I1, I2,
                    soa_1, soa_2,
                    q1, q2, m1, m2, amrex::Real(-1.0), amrex::Real(-1.0),
                    dt, m_CoulombLog, dV, engine );
        }

private:
    amrex::Real m_CoulombLog;
};

#endif // PAIRWISE_COULOMB_COLLISION_FUNC_H_
