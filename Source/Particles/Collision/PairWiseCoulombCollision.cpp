/* Copyright 2020 Yinjian Zhao, David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PairWiseCoulombCollision.H"
#include "ShuffleFisherYates.H"
#include "ElasticCollisionPerez.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

using namespace amrex::literals;

PairWiseCoulombCollision::PairWiseCoulombCollision (std::string const collision_name)
    : CollisionBase(collision_name)
{

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() == 2,
                                     "Pair wise Coulomb must have exactly two species.");

    amrex::ParmParse pp(collision_name);

    // default Coulomb log, if < 0, will be computed automatically
    m_CoulombLog = -1.0_rt;
    queryWithParser(pp, "CoulombLog", m_CoulombLog);

    if (m_species_names[0] == m_species_names[1])
        m_isSameSpecies = true;
    else
        m_isSameSpecies = false;

}

void
PairWiseCoulombCollision::doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc)
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

        // Loop over all grids/tiles at this level
#ifdef _OPENMP
        info.SetDynamic(true);
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi = species1.MakeMFIter(lev, info); mfi.isValid(); ++mfi){
            doCoulombCollisionsWithinTile( lev, mfi, species1, species2 );
        }
    }
}


// Define shortcuts for frequently-used type names
using ParticleType = WarpXParticleContainer::ParticleType;
using ParticleTileType = WarpXParticleContainer::ParticleTileType;
using ParticleBins = amrex::DenseBins<ParticleType>;
using index_type = ParticleBins::index_type;

using namespace ParticleUtils;

/** Perform all binary collisions within a tile
 *
 * @param lev AMR level of the tile
 * @param mfi iterator for multifab
 * @param species1/2 pointer to species container
 * @param ndt user input number of time stpes between collisions
 *
 */
void PairWiseCoulombCollision::doCoulombCollisionsWithinTile
    ( int const lev, amrex::MFIter const& mfi,
    WarpXParticleContainer& species_1,
    WarpXParticleContainer& species_2)
{

    int const ndt = m_ndt;
    amrex::Real CoulombLog = m_CoulombLog;

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
        auto& soa_1 = ptile_1.GetStructOfArrays();
        amrex::ParticleReal * const AMREX_RESTRICT ux_1 =
            soa_1.GetRealData(PIdx::ux).data();
        amrex::ParticleReal * const AMREX_RESTRICT uy_1 =
            soa_1.GetRealData(PIdx::uy).data();
        amrex::ParticleReal * const AMREX_RESTRICT uz_1  =
            soa_1.GetRealData(PIdx::uz).data();
        amrex::ParticleReal const * const AMREX_RESTRICT w_1 =
            soa_1.GetRealData(PIdx::w).data();
        index_type* indices_1 = bins_1.permutationPtr();
        index_type const* cell_offsets_1 = bins_1.offsetsPtr();
        amrex::Real q1 = species_1.getCharge();
        amrex::Real m1 = species_1.getMass();

        const amrex::Real dt = WarpX::GetInstance().getdt(lev);
        amrex::Geometry const& geom = WarpX::GetInstance().Geom(lev);
        amrex::Box const& cbx = mfi.tilebox(amrex::IntVect::TheZeroVector()); //Cell-centered box
        const auto lo = lbound(cbx);
        const auto hi = ubound(cbx);
        int nz = hi.y-lo.y+1;
#if defined WARPX_DIM_XZ
        auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined WARPX_DIM_RZ
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
                if ( cell_stop_1 - cell_start_1 >= 2 )
                {
                    // shuffle
                    ShuffleFisherYates(
                        indices_1, cell_start_1, cell_half_1, engine );

#if defined WARPX_DIM_RZ
                    int ri = (i_cell - i_cell%nz) / nz;
                    auto dV = MathConst::pi*(2.0_rt*ri+1.0_rt)*dr*dr*dz;
#else
                    amrex::ignore_unused(nz);
#endif

                    // Call the function in order to perform collisions
                    ElasticCollisionPerez(
                        cell_start_1, cell_half_1,
                        cell_half_1, cell_stop_1,
                        indices_1, indices_1,
                        ux_1, uy_1, uz_1, ux_1, uy_1, uz_1, w_1, w_1,
                        q1, q1, m1, m1, amrex::Real(-1.0), amrex::Real(-1.0),
                        dt*ndt, CoulombLog, dV, engine );
                }
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
        auto& soa_1 = ptile_1.GetStructOfArrays();
        amrex::ParticleReal * const AMREX_RESTRICT ux_1 =
            soa_1.GetRealData(PIdx::ux).data();
        amrex::ParticleReal * const AMREX_RESTRICT uy_1 =
            soa_1.GetRealData(PIdx::uy).data();
        amrex::ParticleReal * const AMREX_RESTRICT uz_1 =
            soa_1.GetRealData(PIdx::uz).data();
        amrex::ParticleReal const * const AMREX_RESTRICT w_1 =
            soa_1.GetRealData(PIdx::w).data();
        index_type* indices_1 = bins_1.permutationPtr();
        index_type const* cell_offsets_1 = bins_1.offsetsPtr();
        amrex::Real q1 = species_1.getCharge();
        amrex::Real m1 = species_1.getMass();
        // - Species 2
        auto& soa_2 = ptile_2.GetStructOfArrays();
        amrex::Real* ux_2  = soa_2.GetRealData(PIdx::ux).data();
        amrex::Real* uy_2  = soa_2.GetRealData(PIdx::uy).data();
        amrex::Real* uz_2  = soa_2.GetRealData(PIdx::uz).data();
        amrex::Real* w_2   = soa_2.GetRealData(PIdx::w).data();
        index_type* indices_2 = bins_2.permutationPtr();
        index_type const* cell_offsets_2 = bins_2.offsetsPtr();
        amrex::Real q2 = species_2.getCharge();
        amrex::Real m2 = species_2.getMass();

        const amrex::Real dt = WarpX::GetInstance().getdt(lev);
        amrex::Geometry const& geom = WarpX::GetInstance().Geom(lev);
        amrex::Box const& cbx = mfi.tilebox(amrex::IntVect::TheZeroVector()); //Cell-centered box
        const auto lo = lbound(cbx);
        const auto hi = ubound(cbx);
        int nz = hi.y-lo.y+1;
#if defined WARPX_DIM_XZ
        auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined WARPX_DIM_RZ
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
                if ( cell_stop_1 - cell_start_1 >= 1 &&
                     cell_stop_2 - cell_start_2 >= 1 )
                {
                    // shuffle
                    ShuffleFisherYates(indices_1, cell_start_1, cell_stop_1, engine);
                    ShuffleFisherYates(indices_2, cell_start_2, cell_stop_2, engine);

#if defined WARPX_DIM_RZ
                    int ri = (i_cell - i_cell%nz) / nz;
                    auto dV = MathConst::pi*(2.0_rt*ri+1.0_rt)*dr*dr*dz;
#else
                    amrex::ignore_unused(nz);
#endif

                    // Call the function in order to perform collisions
                    ElasticCollisionPerez(
                        cell_start_1, cell_stop_1, cell_start_2, cell_stop_2,
                        indices_1, indices_2,
                        ux_1, uy_1, uz_1, ux_2, uy_2, uz_2, w_1, w_2,
                        q1, q2, m1, m2, amrex::Real(-1.0), amrex::Real(-1.0),
                        dt*ndt, CoulombLog, dV, engine );
                }
            }
        );
    } // end if ( m_isSameSpecies)

}
