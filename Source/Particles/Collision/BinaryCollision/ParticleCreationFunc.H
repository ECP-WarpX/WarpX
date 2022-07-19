/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PARTICLE_CREATION_FUNC_H_
#define PARTICLE_CREATION_FUNC_H_

#include "BinaryCollisionUtils.H"

#include "Particles/Collision/BinaryCollision/NuclearFusion/TwoProductFusionInitializeMomentum.H"
#include "Particles/Collision/BinaryCollision/NuclearFusion/ProtonBoronFusionInitializeMomentum.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_DenseBins.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_INT.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * \brief This functor creates particles produced from a binary collision and sets their initial
 * properties (position, momentum, weight).
 */
class ParticleCreationFunc{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the ParticleCreationFunc class.
     */
    ParticleCreationFunc () = default;

    /**
     * \brief Constructor of the ParticleCreationFunc class
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc pointer to the MultiParticleContainer
     */
    ParticleCreationFunc (const std::string collision_name, MultiParticleContainer const * const mypc);

    /**
     * \brief operator() of the ParticleCreationFunc class. It creates new particles from binary
     * collisions.
     * One product particle is created at the position of each parent particle that collided,
     * allowing for exact charge conservation. For example, in the nuclear reaction
     * "proton + boron -> 3 alpha", we actually create 6 new alpha particles, 3 at the position of
     * of the proton and 3 at the position of the boron.
     * This function also sets the initial weight of the produced particles and subtracts it from
     * the parent particles. If the weight of a parent particle becomes 0, then that particle is
     * deleted.
     * Finally, this function sets the initial momentum of the product particles, by calling a
     * function specific to the considered binary collision.
     *
     * @param[in] n_total_pairs how many binary collisions have been performed in this tile
     * @param[in, out] soa_1 struct of array data of the first colliding particle species
     * @param[in, out] soa_2 struct of array data of the second colliding particle species
     * @param[out] tile_products array containing tile data of the product particles.
     * @param[out] particle_ptr_1 pointer to data of the first colliding particle species. Is
     *             needed to set the id of a particle to -1 in order to delete it when its weight
     *             reaches 0.
     * @param[out] particle_ptr_2 pointer to data of the second colliding particle species. Is
     *             needed to set the id of a particle to -1 in order to delete it when its weight
     *             reaches 0.
     * @param[in] m1 mass of the first colliding particle species
     * @param[in] m2 mass of the second colliding particle species
     * @param[in] products_mass array storing the mass of product particles
     * @param[in] p_mask a mask that is 1 if binary collision has resulted in particle creation
     *            event, 0 otherwise.
     * @param[in] products_np array storing the number of existing product particles in that tile
     * @param[in] copy_species1 array of functors used to copy data from the first colliding
     * particle species to the product particles and to default initialize product particle
     * quantities.
     * @param[in] copy_species2 array of functors used to copy data from the second colliding
     * particle species to the product particles and to default initialize product particle
     * quantities.
     * @param[in] p_pair_indices_1 array that stores the indices of the particles of the first
     * colliding species that were used in the binary collisions (i.e. particle with index
     * p_pair_indices_1[i] took part in collision i)
     * @param[in] p_pair_indices_2 array that stores the indices of the particles of the second
     * colliding species that were used in the binary collisions (i.e. particle with index
     * p_pair_indices_2[i] took part in collision i)
     * @param[in] p_pair_reaction_weight array that stores the weight of the binary collisions.
     * This weight is removed from the parent particles and given to the product particles.
     */
    AMREX_INLINE
    amrex::Vector<int> operator() (
                    const index_type& n_total_pairs,
                    const SoaData_type soa_1, const SoaData_type soa_2,
                    ParticleTileType** AMREX_RESTRICT tile_products,
                    ParticleType* particle_ptr_1, ParticleType*  particle_ptr_2,
                    const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
                    const amrex::Vector<amrex::ParticleReal>& products_mass,
                    const index_type* AMREX_RESTRICT p_mask,
                    const amrex::Vector<index_type>& products_np,
                    const SmartCopy* AMREX_RESTRICT copy_species1,
                    const SmartCopy* AMREX_RESTRICT copy_species2,
                    const index_type* AMREX_RESTRICT p_pair_indices_1,
                    const index_type* AMREX_RESTRICT p_pair_indices_2,
                    const amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight
                    ) const
    {
        if (n_total_pairs == 0) return amrex::Vector<int>(m_num_product_species, 0);

        // Compute offset array and allocate memory for the produced species
        amrex::Gpu::DeviceVector<index_type> offsets(n_total_pairs);
        const auto total = amrex::Scan::ExclusiveSum(n_total_pairs, p_mask, offsets.data());
        const index_type* AMREX_RESTRICT p_offsets = offsets.dataPtr();
        amrex::Vector<int> num_added_vec(m_num_product_species);
        for (int i = 0; i < m_num_product_species; i++)
        {
            // How many particles of product species i are created.
            // Factor 2 is here because we currently create one product species at the position of
            // each source particle of the binary collision. E.g., if a binary collision produces
            // one electron, we create two electrons, one at the position of each particle that
            // collided. This allows for exact charge conservation.
            const index_type num_added = total * m_num_products_host[i] * 2;
            num_added_vec[i] = num_added;
            tile_products[i]->resize(products_np[i] + num_added);
        }

        amrex::ParticleReal* AMREX_RESTRICT w1 = soa_1.m_rdata[PIdx::w];
        amrex::ParticleReal* AMREX_RESTRICT w2 = soa_2.m_rdata[PIdx::w];

        // Create necessary GPU vectors, that will be used in the kernel below
        amrex::Vector<SoaData_type> soa_products;
        for (int i = 0; i < m_num_product_species; i++)
        {
            soa_products.push_back(tile_products[i]->getParticleTileData());
        }
#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<SoaData_type> device_soa_products(m_num_product_species);
        amrex::Gpu::DeviceVector<index_type> device_products_np(m_num_product_species);
        amrex::Gpu::DeviceVector<amrex::ParticleReal> device_products_mass(m_num_product_species);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, soa_products.begin(),
                              soa_products.end(),
                              device_soa_products.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, products_np.begin(),
                              products_np.end(),
                              device_products_np.begin());
       amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, products_mass.begin(),
                              products_mass.end(),
                              device_products_mass.begin());
        amrex::Gpu::streamSynchronize();
        SoaData_type* AMREX_RESTRICT soa_products_data = device_soa_products.data();
        const index_type* AMREX_RESTRICT products_np_data = device_products_np.data();
        const amrex::ParticleReal* AMREX_RESTRICT products_mass_data = device_products_mass.data();
#else
        SoaData_type* AMREX_RESTRICT soa_products_data = soa_products.data();
        const index_type* AMREX_RESTRICT products_np_data = products_np.data();
        const amrex::ParticleReal* AMREX_RESTRICT products_mass_data = products_mass.data();
#endif

        const int t_num_product_species = m_num_product_species;
        const int* AMREX_RESTRICT p_num_products_device = m_num_products_device.data();
        const CollisionType t_collision_type = m_collision_type;

        amrex::ParallelForRNG(n_total_pairs,
        [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
        {
            if (p_mask[i])
            {
                for (int j = 0; j < t_num_product_species; j++)
                {
                    for (int k = 0; k < p_num_products_device[j]; k++)
                    {
                        // Factor 2 is here because we create one product species at the position
                        // of each source particle
                        const auto product_index = products_np_data[j] +
                                                   2*(p_offsets[i]*p_num_products_device[j] + k);
                        // Create product particle at position of particle 1
                        copy_species1[j](soa_products_data[j], soa_1, p_pair_indices_1[i],
                                      product_index, engine);
                        // Create another product particle at position of particle 2
                        copy_species2[j](soa_products_data[j], soa_2, p_pair_indices_2[i],
                                      product_index + 1, engine);

                        // Set the weight of the new particles to p_pair_reaction_weight[i]/2
                        soa_products_data[j].m_rdata[PIdx::w][product_index] =
                                                p_pair_reaction_weight[i]/amrex::ParticleReal(2.);
                        soa_products_data[j].m_rdata[PIdx::w][product_index + 1] =
                                                p_pair_reaction_weight[i]/amrex::ParticleReal(2.);
                    }
                }

                // Remove p_pair_reaction_weight[i] from the colliding particles' weights
                amrex::Gpu::Atomic::AddNoRet(&w1[p_pair_indices_1[i]],
                                                -p_pair_reaction_weight[i]);
                amrex::Gpu::Atomic::AddNoRet(&w2[p_pair_indices_2[i]],
                                                -p_pair_reaction_weight[i]);

                // If the colliding particle weight decreases to zero, remove particle by
                // setting its id to -1
                constexpr amrex::Long minus_one_long = -1;
                if (w1[p_pair_indices_1[i]] <= amrex::ParticleReal(0.))
                {
                    particle_ptr_1[p_pair_indices_1[i]].atomicSetID(minus_one_long);
                }
                if (w2[p_pair_indices_2[i]] <= amrex::ParticleReal(0.))
                {
                    particle_ptr_2[p_pair_indices_2[i]].atomicSetID(minus_one_long);
                }

                // Initialize the product particles' momentum, using a function depending on the
                // specific collision type
                if (t_collision_type == CollisionType::ProtonBoronFusion)
                {
                    const index_type product_start_index = products_np_data[0] + 2*p_offsets[i]*
                                                           p_num_products_device[0];
                    ProtonBoronFusionInitializeMomentum(soa_1, soa_2, soa_products_data[0],
                                                        p_pair_indices_1[i], p_pair_indices_2[i],
                                                        product_start_index, m1, m2, engine);
                }
                else if (t_collision_type == CollisionType::DeuteriumTritiumFusion)
                {
                    amrex::ParticleReal fusion_energy = 0.0_prt;
                    if (t_collision_type == CollisionType::DeuteriumTritiumFusion) {
                        fusion_energy = 17.5893e6_prt * PhysConst::q_e; // 17.6 MeV
                    }
                    TwoProductFusionInitializeMomentum(soa_1, soa_2,
                        soa_products_data[0], soa_products_data[1],
                        p_pair_indices_1[i], p_pair_indices_2[i],
                        products_np_data[0] + 2*p_offsets[i]*p_num_products_device[0],
                        products_np_data[1] + 2*p_offsets[i]*p_num_products_device[1],
                        m1, m2, products_mass_data[0], products_mass_data[1], fusion_energy, engine);
                }

            }
        });

        amrex::Gpu::synchronize();

        return num_added_vec;
    }

private:
    // How many different type of species the collision produces
    int m_num_product_species;
    // Vectors of size m_num_product_species storing how many particles of a given species are
    // produced by a collision event. These vectors are duplicated (one version for host and one
    // for device) which is necessary with GPUs but redundant on CPU.
    amrex::Gpu::DeviceVector<int> m_num_products_device;
    amrex::Gpu::HostVector<int> m_num_products_host;
    CollisionType m_collision_type;
};


/**
 * \brief This class does nothing and is used as second template parameter for binary collisions
 * that do not create particles.
 */
class NoParticleCreationFunc{
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    NoParticleCreationFunc () = default;

    NoParticleCreationFunc (const std::string /*collision_name*/,
                            MultiParticleContainer const * const /*mypc*/) {}

    AMREX_INLINE
    amrex::Vector<int> operator() (
                    const index_type& /*n_total_pairs*/,
                    const SoaData_type /*soa_1*/, const SoaData_type /*soa_2*/,
                    ParticleTileType** /*tile_products*/,
                    ParticleType* /*particle_ptr_1*/, ParticleType* /*particle_ptr_2*/,
                    const amrex::ParticleReal& /*m1*/, const amrex::ParticleReal& /*m2*/,
                    const amrex::Vector<amrex::ParticleReal>& /*products_mass*/,
                    const index_type* /*p_mask*/, const amrex::Vector<index_type>& /*products_np*/,
                    const SmartCopy* /*copy_species1*/, const SmartCopy* /*copy_species2*/,
                    const index_type* /*p_pair_indices_1*/, const index_type* /*p_pair_indices_2*/,
                    const amrex::ParticleReal* /*p_pair_reaction_weight*/
                    ) const
    {
        return amrex::Vector<int>();
    }
};

#endif // PARTICLE_CREATION_FUNC_H_
