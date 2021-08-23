/* Copyright 2020-2021 Yinjian Zhao, David Grote, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_

#include "Particles/Collision/CollisionBase.H"
#include "Particles/Collision/PairWiseCoulombCollisionFunc.H"
#include "Particles/Collision/ShuffleFisherYates.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "WarpX.H"

#include "Particles/MultiParticleContainer_fwd.H"
#include "Particles/WarpXParticleContainer_fwd.H"

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_BLassert.H>
#include <AMReX_Config.H>
#include <AMReX_DenseBins.H>
#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <cmath>
#include <string>

/**
 * \brief This class performs generic binary collisions.
 *
 * \tparam FunctorType the type of the specific binary collision functor that acts on a single cell
 *
 */
template <typename FunctorType>
class BinaryCollision final
    : public CollisionBase
{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;

public:
    /**
     * \brief Constructor of the BinaryCollision class.
     *
     * @param collision_name the name of the collision
     *
     */
    BinaryCollision (std::string collision_name)
        : CollisionBase(collision_name)
    {
        if(m_species_names.size() != 2)
            amrex::Abort("Binary collision " + collision_name + " must have exactly two species.");

        m_binary_collision_functor = FunctorType(collision_name);

        if (m_species_names[0] == m_species_names[1])
            m_isSameSpecies = true;
        else
            m_isSameSpecies = false;
    }

    virtual ~BinaryCollision () = default;

    /** Perform the collisions
     *
     * @param lev AMR level of the tile
     * @param cur_time Current time
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc) override
    {
        const amrex::Real dt = WarpX::GetInstance().getdt(0);
        if ( int(std::floor(cur_time/dt)) % m_ndt != 0 ) return;

        auto& species1 = mypc->GetParticleContainerFromName(m_species_names[0]);
        auto& species2 = mypc->GetParticleContainerFromName(m_species_names[1]);

        // Enable tiling
        amrex::MFItInfo info;
        if (amrex::Gpu::notInLaunchRegion()) info.EnableTiling(species1.tile_size);

        // Loop over refinement levels
        for (int lev = 0; lev <= species1.finestLevel(); ++lev){

        amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

            // Loop over all grids/tiles at this level
#ifdef AMREX_USE_OMP
            info.SetDynamic(true);
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi = species1.MakeMFIter(lev, info); mfi.isValid(); ++mfi){
                if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
                {
                    amrex::Gpu::synchronize();
                }
                amrex::Real wt = amrex::second();

                doCollisionsWithinTile( lev, mfi, species1, species2 );

                if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
                {
                    amrex::Gpu::synchronize();
                    wt = amrex::second() - wt;
                    amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
                }
            }
        }
    }

    /** Perform all binary collisions within a tile
     *
     * @param mfi iterator for multifab
     * @param species1/2 pointer to species container
     *
     */
    void doCollisionsWithinTile (
        int const lev, amrex::MFIter const& mfi,
        WarpXParticleContainer& species_1,
        WarpXParticleContainer& species_2)
    {
        using namespace ParticleUtils;
        using namespace amrex::literals;

        int const ndt = m_ndt;
        FunctorType binary_collision_functor = m_binary_collision_functor;

        if ( m_isSameSpecies ) // species_1 == species_2
        {
            // Extract particles in the tile that `mfi` points to
            ParticleTileType& ptile_1 = species_1.ParticlesAt(lev, mfi);

            // Find the particles that are in each cell of this tile
            ParticleBins bins_1 = findParticlesInEachCell( lev, mfi, ptile_1 );

            // Loop over cells, and collide the particles in each cell

            // Extract low-level data
            int const n_cells = bins_1.numBins();
            // - Species 1
            const auto soa_1 = ptile_1.getParticleTileData();
            index_type* indices_1 = bins_1.permutationPtr();
            index_type const* cell_offsets_1 = bins_1.offsetsPtr();
            amrex::Real q1 = species_1.getCharge();
            amrex::Real m1 = species_1.getMass();
            constexpr int getpos_offset = 0;
            auto get_position_1  = GetParticlePosition(ptile_1, getpos_offset);

            const amrex::Real dt = WarpX::GetInstance().getdt(lev);
            amrex::Geometry const& geom = WarpX::GetInstance().Geom(lev);
#if defined WARPX_DIM_XZ
            auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined WARPX_DIM_RZ
            amrex::Box const& cbx = mfi.tilebox(amrex::IntVect::TheZeroVector()); //Cell-centered box
            const auto lo = lbound(cbx);
            const auto hi = ubound(cbx);
            int const nz = hi.y-lo.y+1;
            auto dr = geom.CellSize(0);
            auto dz = geom.CellSize(1);
#elif (AMREX_SPACEDIM == 3)
            auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

            // Loop over cells
            amrex::ParallelForRNG( n_cells,
                [=] AMREX_GPU_DEVICE (int i_cell, amrex::RandomEngine const& engine) noexcept
                {
                    // The particles from species1 that are in the cell `i_cell` are
                    // given by the `indices_1[cell_start_1:cell_stop_1]`
                    index_type const cell_start_1 = cell_offsets_1[i_cell];
                    index_type const cell_stop_1  = cell_offsets_1[i_cell+1];
                    index_type const cell_half_1 = (cell_start_1+cell_stop_1)/2;

                    // Do not collide if there is only one particle in the cell
                    if ( cell_stop_1 - cell_start_1 <= 1 ) return;

                    // shuffle
                    ShuffleFisherYates(
                        indices_1, cell_start_1, cell_half_1, engine );
#if defined WARPX_DIM_RZ
                    int ri = (i_cell - i_cell%nz) / nz;
                    auto dV = MathConst::pi*(2.0_rt*ri+1.0_rt)*dr*dr*dz;
#endif
                    // Call the function in order to perform collisions
                    binary_collision_functor(
                        cell_start_1, cell_half_1,
                        cell_half_1, cell_stop_1,
                        indices_1, indices_1,
                        soa_1, soa_1, get_position_1, get_position_1,
                        q1, q1, m1, m1, dt*ndt, dV, engine );
                }
            );
        }
        else // species_1 != species_2
        {
            // Extract particles in the tile that `mfi` points to
            ParticleTileType& ptile_1 = species_1.ParticlesAt(lev, mfi);
            ParticleTileType& ptile_2 = species_2.ParticlesAt(lev, mfi);

            // Find the particles that are in each cell of this tile
            ParticleBins bins_1 = findParticlesInEachCell( lev, mfi, ptile_1 );
            ParticleBins bins_2 = findParticlesInEachCell( lev, mfi, ptile_2 );

            // Loop over cells, and collide the particles in each cell

            // Extract low-level data
            int const n_cells = bins_1.numBins();
            // - Species 1
            const auto soa_1 = ptile_1.getParticleTileData();
            index_type* indices_1 = bins_1.permutationPtr();
            index_type const* cell_offsets_1 = bins_1.offsetsPtr();
            amrex::Real q1 = species_1.getCharge();
            amrex::Real m1 = species_1.getMass();
            constexpr int getpos_offset = 0;
            auto get_position_1  = GetParticlePosition(ptile_1, getpos_offset);
            // - Species 2
            const auto soa_2 = ptile_2.getParticleTileData();
            index_type* indices_2 = bins_2.permutationPtr();
            index_type const* cell_offsets_2 = bins_2.offsetsPtr();
            amrex::Real q2 = species_2.getCharge();
            amrex::Real m2 = species_2.getMass();
            auto get_position_2  = GetParticlePosition(ptile_2, getpos_offset);

            const amrex::Real dt = WarpX::GetInstance().getdt(lev);
            amrex::Geometry const& geom = WarpX::GetInstance().Geom(lev);
#if defined WARPX_DIM_XZ
            auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined WARPX_DIM_RZ
            amrex::Box const& cbx = mfi.tilebox(amrex::IntVect::TheZeroVector()); //Cell-centered box
            const auto lo = lbound(cbx);
            const auto hi = ubound(cbx);
            int nz = hi.y-lo.y+1;
            auto dr = geom.CellSize(0);
            auto dz = geom.CellSize(1);
#elif (AMREX_SPACEDIM == 3)
            auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

            // Loop over cells
            amrex::ParallelForRNG( n_cells,
                [=] AMREX_GPU_DEVICE (int i_cell, amrex::RandomEngine const& engine) noexcept
                {
                    // The particles from species1 that are in the cell `i_cell` are
                    // given by the `indices_1[cell_start_1:cell_stop_1]`
                    index_type const cell_start_1 = cell_offsets_1[i_cell];
                    index_type const cell_stop_1  = cell_offsets_1[i_cell+1];
                    // Same for species 2
                    index_type const cell_start_2 = cell_offsets_2[i_cell];
                    index_type const cell_stop_2  = cell_offsets_2[i_cell+1];

                    // ux from species1 can be accessed like this:
                    // ux_1[ indices_1[i] ], where i is between
                    // cell_start_1 (inclusive) and cell_start_2 (exclusive)

                    // Do not collide if one species is missing in the cell
                    if ( cell_stop_1 - cell_start_1 < 1 ||
                        cell_stop_2 - cell_start_2 < 1 ) return;

                    // shuffle
                    ShuffleFisherYates(indices_1, cell_start_1, cell_stop_1, engine);
                    ShuffleFisherYates(indices_2, cell_start_2, cell_stop_2, engine);
#if defined WARPX_DIM_RZ
                    int ri = (i_cell - i_cell%nz) / nz;
                    auto dV = MathConst::pi*(2.0_rt*ri+1.0_rt)*dr*dr*dz;
#endif
                    // Call the function in order to perform collisions
                    binary_collision_functor(
                        cell_start_1, cell_stop_1, cell_start_2, cell_stop_2,
                        indices_1, indices_2,
                        soa_1, soa_2, get_position_1, get_position_2,
                        q1, q2, m1, m2, dt*ndt, dV, engine );
                }
            );
        } // end if ( m_isSameSpecies)

    }

private:

    bool m_isSameSpecies;
    // functor that performs collisions within a cell
    FunctorType m_binary_collision_functor;

};

#endif // WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_
