/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#include "CollisionFilterFunc.H"
#include "DSMC.H"
#include "SplitAndScatterFunc.H"


DSMC::DSMC (const std::string collision_name)
    : CollisionBase(collision_name)
{
    using namespace amrex::literals;
    amrex::ParmParse pp_collision_name(collision_name);

#if defined WARPX_DIM_RZ
    amrex::Abort("DSMC collisions are only implemented for Cartesian coordinates.");
#endif

    if(m_species_names.size() != 2)
    {
        amrex::Abort("DSMC collision " + collision_name + " must have exactly two species.");
    }

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
}

void
DSMC::doCollisions (amrex::Real /*cur_time*/, amrex::Real dt, MultiParticleContainer* mypc)
{
    WARPX_PROFILE("DSMC::doCollisions()");

    auto& species1 = mypc->GetParticleContainerFromName(m_species_names[0]);
    auto& species2 = mypc->GetParticleContainerFromName(m_species_names[1]);

    // SmartCopy objects are created that will facilitate the particle splitting
    // operation involved in DSMC collisions between particles with arbitrary
    // weights.
    SmartCopyFactory copy_factory_species1(species1, species1);
    SmartCopyFactory copy_factory_species2(species2, species2);
    auto copy_species1 = copy_factory_species1.getSmartCopy();
    auto copy_species2 = copy_factory_species2.getSmartCopy();

    species1.defineAllParticleTiles();
    species2.defineAllParticleTiles();

    // Enable tiling
    amrex::MFItInfo info;
    if (amrex::Gpu::notInLaunchRegion()) { info.EnableTiling(WarpXParticleContainer::tile_size); }

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
            auto wt = static_cast<amrex::Real>(amrex::second());

            doCollisionsWithinTile(dt, lev, mfi, species1, species2,
                                   copy_species1, copy_species2);

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = static_cast<amrex::Real>(amrex::second()) - wt;
                amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
            }
        }

        // Call redistribute to remove particles with negative ids
        species1.Redistribute(lev, lev, 0, true, true);
        species2.Redistribute(lev, lev, 0, true, true);
    }
}

void
DSMC::doCollisionsWithinTile(
    amrex::Real dt, int const lev, amrex::MFIter const& mfi,
    WarpXParticleContainer& species_1,
    WarpXParticleContainer& species_2,
    SmartCopy& copy_species1,
    SmartCopy& copy_species2)
{
    using namespace ParticleUtils;
    using namespace amrex::literals;

    // get collision processes
    auto *scattering_processes = m_scattering_processes_exe.data();
    int const process_count   = static_cast<int>(m_scattering_processes_exe.size());

    // Extract particles in the tile that `mfi` points to
    ParticleTileType& ptile_1 = species_1.ParticlesAt(lev, mfi);
    ParticleTileType& ptile_2 = species_2.ParticlesAt(lev, mfi);

    // Find the particles that are in each cell of this tile
    ParticleBins bins_1 = findParticlesInEachCell( lev, mfi, ptile_1 );
    ParticleBins bins_2 = findParticlesInEachCell( lev, mfi, ptile_2 );

    // Extract low-level data
    int const n_cells = static_cast<int>(bins_1.numBins());

    // - Species 1
    index_type* AMREX_RESTRICT indices_1 = bins_1.permutationPtr();
    index_type const* AMREX_RESTRICT cell_offsets_1 = bins_1.offsetsPtr();
    amrex::ParticleReal m1 = species_1.getMass();
    const auto ptd_1 = ptile_1.getParticleTileData();

    // - Species 2
    index_type* AMREX_RESTRICT indices_2 = bins_2.permutationPtr();
    index_type const* AMREX_RESTRICT cell_offsets_2 = bins_2.offsetsPtr();
    amrex::ParticleReal m2 = species_2.getMass();
    const auto ptd_2 = ptile_2.getParticleTileData();

    amrex::Geometry const& geom = WarpX::GetInstance().Geom(lev);
#if defined WARPX_DIM_1D_Z
    auto dV = geom.CellSize(0);
#elif defined WARPX_DIM_XZ
    auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined WARPX_DIM_RZ
    amrex::Box const& cbx = mfi.tilebox(amrex::IntVect::TheZeroVector()); //Cell-centered box
    const auto lo = lbound(cbx);
    const auto hi = ubound(cbx);
    int nz = hi.y-lo.y+1;
    auto dr = geom.CellSize(0);
    auto dz = geom.CellSize(1);
#elif defined(WARPX_DIM_3D)
    auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

    // In the following we set up the "mask" used for creating new particles
    // (from splitting events). There is a mask index for every collision
    // pair. Below we find the size of the mask based on the greater of the
    // number of each species' particle in each cell.
    amrex::Gpu::DeviceVector<index_type> n_pairs_in_each_cell(n_cells);
    index_type* AMREX_RESTRICT p_n_pairs_in_each_cell = n_pairs_in_each_cell.dataPtr();

    // Compute how many pairs in each cell and store in n_pairs_in_each_cell array
    // For different species, the number of pairs in a cell is the number of particles of
    // the species that has the most particles in that cell
    amrex::ParallelFor( n_cells,
        [=] AMREX_GPU_DEVICE (int i_cell) noexcept
        {
            const auto n_part_in_cell_1 = cell_offsets_1[i_cell+1] - cell_offsets_1[i_cell];
            const auto n_part_in_cell_2 = cell_offsets_2[i_cell+1] - cell_offsets_2[i_cell];
            // Particular case: no pair if a species has no particle in that cell
            if (n_part_in_cell_1 == 0 || n_part_in_cell_2 == 0)
            {
                p_n_pairs_in_each_cell[i_cell] = 0;
            }
            else
            {
                p_n_pairs_in_each_cell[i_cell] = amrex::max(n_part_in_cell_1,n_part_in_cell_2);
            }
        }
    );

    // Start indices of the pairs in a cell. Will be used for particle creation
    amrex::Gpu::DeviceVector<index_type> pair_offsets(n_cells);
    const index_type n_total_pairs = (n_cells == 0) ? 0:
                                        amrex::Scan::ExclusiveSum(n_cells,
                                            p_n_pairs_in_each_cell, pair_offsets.data());
    index_type* AMREX_RESTRICT p_pair_offsets = pair_offsets.dataPtr();

    // Now we create the mask. In the DSMC scheme the mask will serve two
    // purposes, 1) record whether a given pair should collide and 2) record
    // the scattering process that should occur.
    amrex::Gpu::DeviceVector<index_type> mask(n_total_pairs);
    index_type* AMREX_RESTRICT p_mask = mask.dataPtr();

    // Will be filled with the index of the first particle of a given pair
    amrex::Gpu::DeviceVector<index_type> pair_indices_1(n_total_pairs);
    index_type* AMREX_RESTRICT p_pair_indices_1 = pair_indices_1.dataPtr();
    // Will be filled with the index of the second particle of a given pair
    amrex::Gpu::DeviceVector<index_type> pair_indices_2(n_total_pairs);
    index_type* AMREX_RESTRICT p_pair_indices_2 = pair_indices_2.dataPtr();

    // How much weight should be given to the produced particle - based on the
    // weight of the collision partner which is not split
    amrex::Gpu::DeviceVector<amrex::ParticleReal> pair_reaction_weight(n_total_pairs);
    amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight =
                                                        pair_reaction_weight.dataPtr();

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
            // Same but for the pairs
            index_type const cell_start_pair = p_pair_offsets[i_cell];

            // ux from species1 can be accessed like this:
            // ux_1[ indices_1[i] ], where i is between
            // cell_start_1 (inclusive) and cell_start_2 (exclusive)

            // Do not collide if one species is missing in the cell
            if ( cell_stop_1 - cell_start_1 < 1 ||
                cell_stop_2 - cell_start_2 < 1 ) { return; }

            // shuffle
            ShuffleFisherYates(indices_1, cell_start_1, cell_stop_1, engine);
            ShuffleFisherYates(indices_2, cell_start_2, cell_stop_2, engine);
#if defined WARPX_DIM_RZ
            int ri = (i_cell - i_cell%nz) / nz;
            auto dV = MathConst::pi*(2.0_prt*ri+1.0_prt)*dr*dr*dz;
#endif
            // Call the function in order to perform collisions
            // If there are product species, p_mask, p_pair_indices_1/2, and
            // p_pair_reaction_weight are filled here
            CollisionFilter(
                cell_start_1, cell_stop_1, cell_start_2, cell_stop_2,
                indices_1, indices_2,
                ptd_1, ptd_2,
                m1, m2, dt, dV,
                cell_start_pair, p_mask, p_pair_indices_1, p_pair_indices_2,
                p_pair_reaction_weight,
                process_count, scattering_processes, engine );
        }
    );

    const auto num_p_tile1 = ptile_1.numParticles();
    const auto num_p_tile2 = ptile_2.numParticles();

    // Create the new product particles and define their initial values
    // num_added: how many particles of each product species have been created
    const int num_added = splitScatteringParticles(
        n_total_pairs,
        ptile_1, ptile_2,
        p_mask,
        copy_species1, copy_species2,
        m1, m2,
        p_pair_indices_1, p_pair_indices_2,
        p_pair_reaction_weight);

    if (num_added > 0) {
        setNewParticleIDs(ptile_1, num_p_tile1, num_added);
        setNewParticleIDs(ptile_2, num_p_tile2, num_added);
    }
}
