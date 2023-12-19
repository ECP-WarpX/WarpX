/* Copyright 2023 The WarpX Community
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
 * \brief This class performs DSMC (direct simulation Monte Carlo) collisions
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
     * \brief Constructor of the DSMC class
     *
     * @param[in] collision_name the name of the collision
     */
    DSMC (std::string collision_name);

    /** Perform the collisions
     *
     * @param cur_time Current time
     * @param dt Time step size
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real /*cur_time*/, amrex::Real dt, MultiParticleContainer* mypc) override;

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
        amrex::Real dt, int lev, amrex::MFIter const& mfi,
        WarpXParticleContainer& species_1,
        WarpXParticleContainer& species_2,
        SmartCopy& copy_species1,
        SmartCopy& copy_species2 );

private:
    amrex::Vector<ScatteringProcess> m_scattering_processes;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor> m_scattering_processes_exe;
};

#endif // DSMC_H_
