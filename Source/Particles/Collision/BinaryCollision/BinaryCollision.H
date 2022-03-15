/* Copyright 2020-2021 Yinjian Zhao, David Grote, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_

#include "Particles/Collision/BinaryCollision/Coulomb/PairWiseCoulombCollisionFunc.H"
#include "Particles/Collision/BinaryCollision/NuclearFusion/NuclearFusionFunc.H"
#include "Particles/Collision/BinaryCollision/ParticleCreationFunc.H"
#include "Particles/Collision/BinaryCollision/ShuffleFisherYates.H"
#include "Particles/Collision/CollisionBase.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
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
#include <AMReX_GpuContainers.H>
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
#include <AMReX_Scan.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <cmath>
#include <string>

/**
 * \brief This class performs generic binary collisions.
 *
 * \tparam CollisionFunctorType the type of the specific binary collision functor that acts on a
 *         single cell
 * \tparam CopyTransformFunctorType the type of the second functor used in the case of
 *         particle creation
 *
 */
template <typename CollisionFunctorType,
          typename CopyTransformFunctorType = NoParticleCreationFunc>
class BinaryCollision final
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
     * \brief Constructor of the BinaryCollision class.
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc Container of species involved
     *
     */
    BinaryCollision (std::string collision_name, MultiParticleContainer const * const mypc)
        : CollisionBase(collision_name)
    {
        if(m_species_names.size() != 2)
            amrex::Abort("Binary collision " + collision_name + " must have exactly two species.");

        if (m_species_names[0] == m_species_names[1])
            m_isSameSpecies = true;
        else
            m_isSameSpecies = false;

        m_binary_collision_functor = CollisionFunctorType(collision_name, mypc, m_isSameSpecies);

        amrex::ParmParse pp_collision_name(collision_name);
        pp_collision_name.queryarr("product_species", m_product_species);
        m_have_product_species = m_product_species.size() > 0;
        m_copy_transform_functor = CopyTransformFunctorType(collision_name, mypc);
    }

    virtual ~BinaryCollision () = default;

    /** Perform the collisions
     *
     * @param cur_time Current time
     * @param dt Time step size
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real /*cur_time*/, amrex::Real dt, MultiParticleContainer* mypc) override
    {

        auto& species1 = mypc->GetParticleContainerFromName(m_species_names[0]);
        auto& species2 = mypc->GetParticleContainerFromName(m_species_names[1]);

        // In case of particle creation, create the necessary vectors
        const int n_product_species = m_product_species.size();
        amrex::Vector<WarpXParticleContainer*> product_species_vector;
        amrex::Vector<SmartCopyFactory> copy_factory_species1;
        amrex::Vector<SmartCopyFactory> copy_factory_species2;
        amrex::Vector<SmartCopy> copy_species1;
        amrex::Vector<SmartCopy> copy_species2;
        for (int i = 0; i < n_product_species; i++)
        {
            auto& product = mypc->GetParticleContainerFromName(m_product_species[i]);
            product.defineAllParticleTiles();
            product_species_vector.push_back(&product);
            // Although the copy factories are not explicitly reused past this point, we need to
            // store them in vectors so that the data that they own, which is used by the smart
            // copy functors, does not go out of scope at the end of this for loop.
            copy_factory_species1.push_back(SmartCopyFactory(species1, product));
            copy_factory_species2.push_back(SmartCopyFactory(species2, product));
            copy_species1.push_back(copy_factory_species1[i].getSmartCopy());
            copy_species2.push_back(copy_factory_species2[i].getSmartCopy());
        }
#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<SmartCopy> device_copy_species1(n_product_species);
        amrex::Gpu::DeviceVector<SmartCopy> device_copy_species2(n_product_species);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, copy_species1.begin(),
                              copy_species1.end(), device_copy_species1.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, copy_species2.begin(),
                              copy_species2.end(), device_copy_species2.begin());
        amrex::Gpu::streamSynchronize();
        auto copy_species1_data = device_copy_species1.data();
        auto copy_species2_data = device_copy_species2.data();
#else
        auto copy_species1_data = copy_species1.data();
        auto copy_species2_data = copy_species2.data();
#endif
        if (m_have_product_species){
            species1.defineAllParticleTiles();
            species2.defineAllParticleTiles();
        }

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

                doCollisionsWithinTile( dt, lev, mfi, species1, species2, product_species_vector,
                                         copy_species1_data, copy_species2_data);

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
     * \param[in] lev the mesh-refinement level
     * \param[in] mfi iterator for multifab
     * \param species_1 first species container
     * \param species_2 second species container
     * \param product_species_vector vector of pointers to product species containers
     * \param copy_species1 vector of SmartCopy functors used to copy species 1 to product species
     * \param copy_species2 vector of SmartCopy functors used to copy species 2 to product species
     *
     */
    void doCollisionsWithinTile (
        amrex::Real dt, int const lev, amrex::MFIter const& mfi,
        WarpXParticleContainer& species_1,
        WarpXParticleContainer& species_2,
        amrex::Vector<WarpXParticleContainer*> product_species_vector,
        SmartCopy* copy_species1, SmartCopy* copy_species2)
    {
        using namespace ParticleUtils;
        using namespace amrex::literals;

        CollisionFunctorType binary_collision_functor = m_binary_collision_functor;
        const bool have_product_species = m_have_product_species;

        // Store product species data in vectors
        const int n_product_species = m_product_species.size();
        amrex::Vector<ParticleTileType*> tile_products;
        amrex::Vector<GetParticlePosition> get_position_products;
        amrex::Vector<index_type> products_np;
        constexpr int getpos_offset = 0;
        for (int i = 0; i < n_product_species; i++)
        {
            ParticleTileType& ptile_product = product_species_vector[i]->ParticlesAt(lev, mfi);
            tile_products.push_back(&ptile_product);
            get_position_products.push_back(GetParticlePosition(ptile_product,
                                                                      getpos_offset));
            products_np.push_back(ptile_product.numParticles());
        }
        auto tile_products_data = tile_products.data();

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
            index_type* AMREX_RESTRICT indices_1 = bins_1.permutationPtr();
            index_type const* AMREX_RESTRICT cell_offsets_1 = bins_1.offsetsPtr();
            amrex::ParticleReal q1 = species_1.getCharge();
            amrex::ParticleReal m1 = species_1.getMass();
            auto get_position_1  = GetParticlePosition(ptile_1, getpos_offset);
            // Needed to access the particle id
            ParticleType * AMREX_RESTRICT
                                 particle_ptr_1 = ptile_1.GetArrayOfStructs()().data();

            amrex::Geometry const& geom = WarpX::GetInstance().Geom(lev);
#if defined WARPX_DIM_1D_Z
            auto dV = geom.CellSize(0);
#elif defined WARPX_DIM_XZ
            auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined WARPX_DIM_RZ
            amrex::Box const& cbx = mfi.tilebox(amrex::IntVect::TheZeroVector()); //Cell-centered box
            const auto lo = lbound(cbx);
            const auto hi = ubound(cbx);
            int const nz = hi.y-lo.y+1;
            auto dr = geom.CellSize(0);
            auto dz = geom.CellSize(1);
#elif defined(WARPX_DIM_3D)
            auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif


            /*
              The following calculations are only required when creating product particles
            */
            const int n_cells_products = have_product_species ? n_cells: 0;
            amrex::Gpu::DeviceVector<index_type> n_pairs_in_each_cell(n_cells_products);
            index_type* AMREX_RESTRICT p_n_pairs_in_each_cell = n_pairs_in_each_cell.dataPtr();

            // Compute how many pairs in each cell and store in n_pairs_in_each_cell array
            // For a single species, the number of pair in a cell is half the number of particles
            // in that cell, rounded up to the next higher integer.
            amrex::ParallelFor( n_cells_products,
                [=] AMREX_GPU_DEVICE (int i_cell) noexcept
                {
                    const auto n_part_in_cell = cell_offsets_1[i_cell+1] - cell_offsets_1[i_cell];
                    // Particular case: if there's only 1 particle in a cell, then there's no pair
                    p_n_pairs_in_each_cell[i_cell] = (n_part_in_cell == 1)? 0: (n_part_in_cell+1)/2;
                }
            );

            // Start indices of the pairs in a cell. Will be used for particle creation.
            amrex::Gpu::DeviceVector<index_type> pair_offsets(n_cells_products);
            const index_type n_total_pairs = (n_cells_products == 0) ? 0:
                                                amrex::Scan::ExclusiveSum(n_cells_products,
                                                    p_n_pairs_in_each_cell, pair_offsets.data());
            index_type* AMREX_RESTRICT p_pair_offsets = pair_offsets.dataPtr();

            // mask: equal to 1 if particle creation occurs for a given pair, 0 otherwise
            amrex::Gpu::DeviceVector<index_type> mask(n_total_pairs);
            index_type* AMREX_RESTRICT p_mask = mask.dataPtr();
            // Will be filled with the index of the first particle of a given pair
            amrex::Gpu::DeviceVector<index_type> pair_indices_1(n_total_pairs);
            index_type* AMREX_RESTRICT p_pair_indices_1 = pair_indices_1.dataPtr();
            // Will be filled with the index of the second particle of a given pair
            amrex::Gpu::DeviceVector<index_type> pair_indices_2(n_total_pairs);
            index_type* AMREX_RESTRICT p_pair_indices_2 = pair_indices_2.dataPtr();
            // How much weight should be given to the produced particles (and removed from the
            // reacting particles)
            amrex::Gpu::DeviceVector<amrex::ParticleReal> pair_reaction_weight(n_total_pairs);
            amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight =
                                                               pair_reaction_weight.dataPtr();
            /*
              End of calculations only required when creating product particles
            */


            // Loop over cells
            amrex::ParallelForRNG( n_cells,
                [=] AMREX_GPU_DEVICE (int i_cell, amrex::RandomEngine const& engine) noexcept
                {
                    // The particles from species1 that are in the cell `i_cell` are
                    // given by the `indices_1[cell_start_1:cell_stop_1]`
                    index_type const cell_start_1 = cell_offsets_1[i_cell];
                    index_type const cell_stop_1  = cell_offsets_1[i_cell+1];
                    index_type const cell_half_1 = (cell_start_1+cell_stop_1)/2;

                    // Same but for the pairs
                    index_type const cell_start_pair = have_product_species?
                                                           p_pair_offsets[i_cell] : 0;

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
                    // If there are product species, mask, p_pair_indices_1/2, and
                    // p_pair_reaction_weight are filled here
                    binary_collision_functor(
                        cell_start_1, cell_half_1,
                        cell_half_1, cell_stop_1,
                        indices_1, indices_1,
                        soa_1, soa_1, get_position_1, get_position_1,
                        q1, q1, m1, m1, dt, dV,
                        cell_start_pair, p_mask, p_pair_indices_1, p_pair_indices_2,
                        p_pair_reaction_weight, engine );
                }
            );

            // Create the new product particles and define their initial values
            // num_added: how many particles of each product species have been created
            const amrex::Vector<int> num_added = m_copy_transform_functor(n_total_pairs,
                                                    soa_1, soa_1, tile_products_data,
                                                    particle_ptr_1, particle_ptr_1, m1, m1,
                                                    p_mask, products_np,
                                                    copy_species1, copy_species2,
                                                    p_pair_indices_1, p_pair_indices_2,
                                                    p_pair_reaction_weight);

            for (int i = 0; i < n_product_species; i++)
            {
                setNewParticleIDs(*(tile_products_data[i]), products_np[i], num_added[i]);
            }
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
            index_type* AMREX_RESTRICT indices_1 = bins_1.permutationPtr();
            index_type const* AMREX_RESTRICT cell_offsets_1 = bins_1.offsetsPtr();
            amrex::ParticleReal q1 = species_1.getCharge();
            amrex::ParticleReal m1 = species_1.getMass();
            auto get_position_1  = GetParticlePosition(ptile_1, getpos_offset);
            // Needed to access the particle id
            ParticleType * AMREX_RESTRICT
                                 particle_ptr_1 = ptile_1.GetArrayOfStructs()().data();
            // - Species 2
            const auto soa_2 = ptile_2.getParticleTileData();
            index_type* AMREX_RESTRICT indices_2 = bins_2.permutationPtr();
            index_type const* AMREX_RESTRICT cell_offsets_2 = bins_2.offsetsPtr();
            amrex::ParticleReal q2 = species_2.getCharge();
            amrex::ParticleReal m2 = species_2.getMass();
            auto get_position_2  = GetParticlePosition(ptile_2, getpos_offset);
            // Needed to access the particle id
            ParticleType * AMREX_RESTRICT
                                 particle_ptr_2 = ptile_2.GetArrayOfStructs()().data();

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


            /*
              The following calculations are only required when creating product particles
            */
            const int n_cells_products = have_product_species ? n_cells: 0;
            amrex::Gpu::DeviceVector<index_type> n_pairs_in_each_cell(n_cells_products);
            index_type* AMREX_RESTRICT p_n_pairs_in_each_cell = n_pairs_in_each_cell.dataPtr();

            // Compute how many pairs in each cell and store in n_pairs_in_each_cell array
            // For different species, the number of pairs in a cell is the number of particles of
            // the species that has the most particles in that cell
            amrex::ParallelFor( n_cells_products,
                [=] AMREX_GPU_DEVICE (int i_cell) noexcept
                {
                    const auto n_part_in_cell_1 = cell_offsets_1[i_cell+1] - cell_offsets_1[i_cell];
                    const auto n_part_in_cell_2 = cell_offsets_2[i_cell+1] - cell_offsets_2[i_cell];
                    // Particular case: no pair if a species has no particle in that cell
                    if (n_part_in_cell_1 == 0 || n_part_in_cell_2 == 0)
                        p_n_pairs_in_each_cell[i_cell] = 0;
                    else
                        p_n_pairs_in_each_cell[i_cell] =
                                                      amrex::max(n_part_in_cell_1,n_part_in_cell_2);
                }
            );

            // Start indices of the pairs in a cell. Will be used for particle creation
            amrex::Gpu::DeviceVector<index_type> pair_offsets(n_cells_products);
            const index_type n_total_pairs = (n_cells_products == 0) ? 0:
                                                amrex::Scan::ExclusiveSum(n_cells_products,
                                                    p_n_pairs_in_each_cell, pair_offsets.data());
            index_type* AMREX_RESTRICT p_pair_offsets = pair_offsets.dataPtr();

            // mask: equal to 1 if particle creation occurs for a given pair, 0 otherwise
            amrex::Gpu::DeviceVector<index_type> mask(n_total_pairs);
            index_type* AMREX_RESTRICT p_mask = mask.dataPtr();
            // Will be filled with the index of the first particle of a given pair
            amrex::Gpu::DeviceVector<index_type> pair_indices_1(n_total_pairs);
            index_type* AMREX_RESTRICT p_pair_indices_1 = pair_indices_1.dataPtr();
            // Will be filled with the index of the second particle of a given pair
            amrex::Gpu::DeviceVector<index_type> pair_indices_2(n_total_pairs);
            index_type* AMREX_RESTRICT p_pair_indices_2 = pair_indices_2.dataPtr();
            // How much weight should be given to the produced particles (and removed from the
            // reacting particles)
            amrex::Gpu::DeviceVector<amrex::ParticleReal> pair_reaction_weight(n_total_pairs);
            amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight =
                                                               pair_reaction_weight.dataPtr();
            /*
              End of calculations only required when creating product particles
            */


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
                    index_type const cell_start_pair = have_product_species?
                                                           p_pair_offsets[i_cell] : 0;

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
                    // If there are product species, p_mask, p_pair_indices_1/2, and
                    // p_pair_reaction_weight are filled here
                    binary_collision_functor(
                        cell_start_1, cell_stop_1, cell_start_2, cell_stop_2,
                        indices_1, indices_2,
                        soa_1, soa_2, get_position_1, get_position_2,
                        q1, q2, m1, m2, dt, dV,
                        cell_start_pair, p_mask, p_pair_indices_1, p_pair_indices_2,
                        p_pair_reaction_weight, engine );
                }
            );


            // Create the new product particles and define their initial values
            // num_added: how many particles of each product species have been created
            const amrex::Vector<int> num_added = m_copy_transform_functor(n_total_pairs,
                                                    soa_1, soa_2, tile_products_data,
                                                    particle_ptr_1, particle_ptr_2, m1, m2,
                                                    p_mask, products_np,
                                                    copy_species1, copy_species2,
                                                    p_pair_indices_1, p_pair_indices_2,
                                                    p_pair_reaction_weight);

            for (int i = 0; i < n_product_species; i++)
            {
                setNewParticleIDs(*(tile_products_data[i]), products_np[i], num_added[i]);
            }

        } // end if ( m_isSameSpecies)

    }

private:

    bool m_isSameSpecies;
    bool m_have_product_species;
    amrex::Vector<std::string> m_product_species;
    // functor that performs collisions within a cell
    CollisionFunctorType m_binary_collision_functor;
    // functor that creates new particles and initializes their parameters
    CopyTransformFunctorType m_copy_transform_functor;

};

#endif // WARPX_PARTICLES_COLLISION_BINARYCOLLISION_H_
