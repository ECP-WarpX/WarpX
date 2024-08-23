/* Copyright 2023-2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPLIT_AND_SCATTER_FUNC_H_
#define WARPX_SPLIT_AND_SCATTER_FUNC_H_

#include "Particles/Collision/BinaryCollision/BinaryCollisionUtils.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/ParticleUtils.H"

/**
 * \brief This class defines an operator to create product particles from DSMC
 * collisions and sets the particle properties (position, momentum, weight).
 */
class SplitAndScatterFunc
{
    // Define shortcuts for frequently-used type names
    using ParticleType = typename WarpXParticleContainer::ParticleType;
    using ParticleTileType = typename WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = typename ParticleTileType::ParticleTileDataType;
    using ParticleBins = amrex::DenseBins<ParticleTileDataType>;
    using index_type = typename ParticleBins::index_type;
    using SoaData_type = typename WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the SplitAndScatterFunc class.
     */
    SplitAndScatterFunc () = default;

    /**
     * \brief Constructor of the SplitAndScatterFunc class
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc pointer to the MultiParticleContainer
     */
    SplitAndScatterFunc (const std::string& collision_name, MultiParticleContainer const * mypc);

    /**
     * \brief Function that performs the particle scattering and injection due
     * to binary collisions.
     *
     * \return num_added the number of particles added to each species.
     */
    AMREX_INLINE
    amrex::Vector<int> operator() (
        const index_type& n_total_pairs,
        ParticleTileType& ptile1, ParticleTileType& ptile2,
        const amrex::Vector<WarpXParticleContainer*>& pc_products,
        ParticleTileType** AMREX_RESTRICT tile_products,
        const amrex::ParticleReal m1, const amrex::ParticleReal m2,
        const amrex::Vector<amrex::ParticleReal>& /*products_mass*/,
        const index_type* AMREX_RESTRICT mask,
        const amrex::Vector<index_type>& products_np,
        const SmartCopy* AMREX_RESTRICT copy_species1,
        const SmartCopy* AMREX_RESTRICT copy_species2,
        const index_type* AMREX_RESTRICT p_pair_indices_1,
        const index_type* AMREX_RESTRICT p_pair_indices_2,
        const amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight ) const
    {
        using namespace amrex::literals;

        if (n_total_pairs == 0) { return amrex::Vector<int>(m_num_product_species, 0); }

        const amrex::ParticleReal min_splitting_weight = m_min_splitting_weight;

        amrex::Gpu::DeviceVector<index_type> offsets(n_total_pairs);
        index_type* AMREX_RESTRICT offsets_data = offsets.data();
        const index_type* AMREX_RESTRICT p_offsets = offsets.dataPtr();

        // The following is used to calculate the appropriate offsets. Note that
        // a standard cummulative sum is not appropriate since the mask is also
        // used to specify the type of collision and can therefore have values >1
        auto const total = amrex::Scan::PrefixSum<index_type>(n_total_pairs,
            [=] AMREX_GPU_DEVICE (index_type i) -> index_type { return (mask[i] && (p_pair_reaction_weight[i] > min_splitting_weight)) ? 1 : 0; },
            [=] AMREX_GPU_DEVICE (index_type i, index_type s) { offsets_data[i] = s; },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum
        );

        amrex::Vector<int> num_added_vec(m_num_product_species);
        for (int i = 0; i < m_num_product_species; i++)
        {
            // How many particles of product species i are created.
            const index_type num_added = total * m_num_products_host[i];
            num_added_vec[i] = static_cast<int>(num_added);
            tile_products[i]->resize(products_np[i] + num_added);
        }

        const auto soa_1 = ptile1.getParticleTileData();
        const auto soa_2 = ptile2.getParticleTileData();

        // Create necessary GPU vectors, that will be used in the kernel below
        amrex::Vector<SoaData_type> soa_products;
        for (int i = 0; i < m_num_product_species; i++)
        {
            soa_products.push_back(tile_products[i]->getParticleTileData());
        }
#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<SoaData_type> device_soa_products(m_num_product_species);
        amrex::Gpu::DeviceVector<index_type> device_products_np(m_num_product_species);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, soa_products.begin(),
                              soa_products.end(),
                              device_soa_products.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, products_np.begin(),
                              products_np.end(),
                              device_products_np.begin());

        amrex::Gpu::streamSynchronize();
        SoaData_type* AMREX_RESTRICT soa_products_data = device_soa_products.data();
        const index_type* AMREX_RESTRICT products_np_data = device_products_np.data();
#else
        SoaData_type* AMREX_RESTRICT soa_products_data = soa_products.data();
        const index_type* AMREX_RESTRICT products_np_data = products_np.data();
#endif

        const int* AMREX_RESTRICT p_num_products_device = m_num_products_device.data();

        constexpr auto c2 = PhysConst::c * PhysConst::c;

        amrex::ParallelForRNG(n_total_pairs,
        [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
        {
            if (mask[i])
            {
                // for now we ignore the possibility of having actual reaction
                // products - only duplicating (splitting) of the colliding
                // particles is supported. If the collision weight is below the
                // user-set minimum weight, the product particles are immediately
                // re-merged with the colliding particles.

                // Get the colliding velocities
                auto ux1 = soa_1.m_rdata[PIdx::ux][p_pair_indices_1[i]];
                auto uy1 = soa_1.m_rdata[PIdx::uy][p_pair_indices_1[i]];
                auto uz1 = soa_1.m_rdata[PIdx::uz][p_pair_indices_1[i]];
                auto ux2 = soa_2.m_rdata[PIdx::ux][p_pair_indices_2[i]];
                auto uy2 = soa_2.m_rdata[PIdx::uy][p_pair_indices_2[i]];
                auto uz2 = soa_2.m_rdata[PIdx::uz][p_pair_indices_2[i]];

                // for simplicity (for now) we assume non-relativistic particles
                // and simply calculate the center-of-momentum velocity from the
                // rest masses
                auto const uCOM_x = (m1 * ux1 + m2 * ux2) / (m1 + m2);
                auto const uCOM_y = (m1 * uy1 + m2 * uy2) / (m1 + m2);
                auto const uCOM_z = (m1 * uz1 + m2 * uz2) / (m1 + m2);

                // transform to COM frame
                ux1 -= uCOM_x;
                uy1 -= uCOM_y;
                uz1 -= uCOM_z;
                ux2 -= uCOM_x;
                uy2 -= uCOM_y;
                uz2 -= uCOM_z;

                if (mask[i] == int(ScatteringProcessType::ELASTIC)) {
                    // randomly rotate the velocity vector for the first particle
                    ParticleUtils::RandomizeVelocity(
                        ux1, uy1, uz1, std::sqrt(ux1*ux1 + uy1*uy1 + uz1*uz1), engine
                    );
                    // set the second particles velocity so that the total momentum
                    // is zero
                    ux2 = -ux1 * m1 / m2;
                    uy2 = -uy1 * m1 / m2;
                    uz2 = -uz1 * m1 / m2;
                } else if (mask[i] == int(ScatteringProcessType::BACK)) {
                    // reverse the velocity vectors of both particles
                    ux1 *= -1.0_prt;
                    uy1 *= -1.0_prt;
                    uz1 *= -1.0_prt;
                    ux2 *= -1.0_prt;
                    uy2 *= -1.0_prt;
                    uz2 *= -1.0_prt;
                } else if (mask[i] == int(ScatteringProcessType::CHARGE_EXCHANGE)) {
                    if (std::abs(m1 - m2) < 1e-28) {
                        auto const temp_ux = ux1;
                        auto const temp_uy = uy1;
                        auto const temp_uz = uz1;
                        ux1 = ux2;
                        uy1 = uy2;
                        uz1 = uz2;
                        ux2 = temp_ux;
                        uy2 = temp_uy;
                        uz2 = temp_uz;
                    }
                    else {
                        amrex::Abort("Uneven mass charge-exchange not implemented yet.");
                    }
                }
                else {
                    amrex::Abort("Unknown scattering process.");
                }
                // transform back to labframe
                ux1 += uCOM_x;
                uy1 += uCOM_y;
                uz1 += uCOM_z;
                ux2 += uCOM_x;
                uy2 += uCOM_y;
                uz2 += uCOM_z;

                // If the reaction weight is more than the minimum set weight
                // split the colliding weight off the colliding particles
                if (p_pair_reaction_weight[i] > min_splitting_weight) {
                    const auto product1_index = products_np_data[0] +
                                            (p_offsets[i]*p_num_products_device[0] + 0);
                    // Make a copy of the particle from species 1
                    copy_species1[0](soa_products_data[0], soa_1, static_cast<int>(p_pair_indices_1[i]),
                                    static_cast<int>(product1_index), engine);
                    // Set the weight of the new particles to p_pair_reaction_weight[i]
                    soa_products_data[0].m_rdata[PIdx::w][product1_index] = p_pair_reaction_weight[i];

                    const auto product2_index = products_np_data[1] +
                                            (p_offsets[i]*p_num_products_device[1] + 0);
                    // Make a copy of the particle from species 2
                    copy_species2[1](soa_products_data[1], soa_2, static_cast<int>(p_pair_indices_2[i]),
                                    static_cast<int>(product2_index), engine);
                    // Set the weight of the new particles to p_pair_reaction_weight[i]
                    soa_products_data[1].m_rdata[PIdx::w][product2_index] = p_pair_reaction_weight[i];

                    // set the new particles' velocties to the post collision velocities
                    soa_products_data[0].m_rdata[PIdx::ux][product1_index] = ux1;
                    soa_products_data[0].m_rdata[PIdx::uy][product1_index] = uy1;
                    soa_products_data[0].m_rdata[PIdx::uz][product1_index] = uz1;
                    soa_products_data[1].m_rdata[PIdx::ux][product2_index] = ux2;
                    soa_products_data[1].m_rdata[PIdx::uy][product2_index] = uy2;
                    soa_products_data[1].m_rdata[PIdx::uz][product2_index] = uz2;

                }
                else {
                    // Modify the colliding particles' velocities without splitting.
                    // The new velocities are rotated based on the scattered velocities
                    // but conserves energy.

                    const auto w1_frac = p_pair_reaction_weight[i] / w1[p_pair_indices_1[i]];
                    const auto w2_frac = p_pair_reaction_weight[i] / w2[p_pair_indices_2[i]];

                    // Calculate energies of the parts of the colliding pairs
                    // that did not participate in the collisions
                    double Ek1, Ek2;
                    ParticleUtils::getEnergy(
                        soa_1.m_rdata[PIdx::ux][p_pair_indices_1[i]],
                        soa_1.m_rdata[PIdx::uy][p_pair_indices_1[i]],
                        soa_1.m_rdata[PIdx::uz][p_pair_indices_1[i]], m1, Ek1
                    );
                    Ek1 *= (1.0_prt - w1_frac);
                    ParticleUtils::getEnergy(
                        soa_2.m_rdata[PIdx::ux][p_pair_indices_2[i]],
                        soa_2.m_rdata[PIdx::uy][p_pair_indices_2[i]],
                        soa_2.m_rdata[PIdx::uz][p_pair_indices_2[i]], m2, Ek2
                    );
                    Ek2 *= (1.0_prt - w2_frac);

                    // Calculate energy of the product particles
                    double Ek_product1, Ek_product2;
                    ParticleUtils::getEnergy(ux1, uy1, uz1, m1, Ek_product1);
                    Ek_product1 *= w1_frac;
                    ParticleUtils::getEnergy(ux2, uy2, uz2, m2, Ek_product2);
                    Ek_product2 *= w2_frac;

                    // Scatter velocities based on weighted sum
                    soa_1.m_rdata[PIdx::ux][p_pair_indices_1[i]] *= 1.0_prt - w1_frac;
                    soa_1.m_rdata[PIdx::uy][p_pair_indices_1[i]] *= 1.0_prt - w1_frac;
                    soa_1.m_rdata[PIdx::uz][p_pair_indices_1[i]] *= 1.0_prt - w1_frac;
                    soa_1.m_rdata[PIdx::ux][p_pair_indices_1[i]] += w1_frac * ux1;
                    soa_1.m_rdata[PIdx::uy][p_pair_indices_1[i]] += w1_frac * uy1;
                    soa_1.m_rdata[PIdx::uz][p_pair_indices_1[i]] += w1_frac * uz1;

                    soa_2.m_rdata[PIdx::ux][p_pair_indices_2[i]] *= 1.0_prt - w2_frac;
                    soa_2.m_rdata[PIdx::uy][p_pair_indices_2[i]] *= 1.0_prt - w2_frac;
                    soa_2.m_rdata[PIdx::uz][p_pair_indices_2[i]] *= 1.0_prt - w2_frac;
                    soa_2.m_rdata[PIdx::ux][p_pair_indices_2[i]] += w2_frac * ux2;
                    soa_2.m_rdata[PIdx::uy][p_pair_indices_2[i]] += w2_frac * uy2;
                    soa_2.m_rdata[PIdx::uz][p_pair_indices_2[i]] += w2_frac * uz2;

                    // rescale velocity vectors to conserve kinetic energy
                    const amrex::ParticleReal up1_mag = std::sqrt(
                        soa_1.m_rdata[PIdx::ux][p_pair_indices_1[i]]*soa_1.m_rdata[PIdx::ux][p_pair_indices_1[i]]
                        + soa_1.m_rdata[PIdx::uy][p_pair_indices_1[i]]*soa_1.m_rdata[PIdx::uy][p_pair_indices_1[i]]
                        + soa_1.m_rdata[PIdx::uz][p_pair_indices_1[i]]*soa_1.m_rdata[PIdx::uz][p_pair_indices_1[i]]
                    );
                    amrex::ParticleReal scale_fac1 = 1.0_prt;
                    if (up1_mag > 0.0_prt) {
                        const auto Etot1 = (Ek1 + Ek_product1) * PhysConst::q_e;
                        scale_fac1 = static_cast<amrex::ParticleReal>(
                            std::sqrt(Etot1 * (Etot1 + 2.0_prt*m1*c2) / c2) / m1
                        ) / up1_mag;
                    }
                    soa_1.m_rdata[PIdx::ux][p_pair_indices_1[i]] *= scale_fac1;
                    soa_1.m_rdata[PIdx::uy][p_pair_indices_1[i]] *= scale_fac1;
                    soa_1.m_rdata[PIdx::uz][p_pair_indices_1[i]] *= scale_fac1;

                    const amrex::ParticleReal up2_mag = std::sqrt(
                        soa_2.m_rdata[PIdx::ux][p_pair_indices_2[i]]*soa_2.m_rdata[PIdx::ux][p_pair_indices_2[i]]
                        + soa_2.m_rdata[PIdx::uy][p_pair_indices_2[i]]*soa_2.m_rdata[PIdx::uy][p_pair_indices_2[i]]
                        + soa_2.m_rdata[PIdx::uz][p_pair_indices_2[i]]*soa_2.m_rdata[PIdx::uz][p_pair_indices_2[i]]
                    );
                    amrex::ParticleReal scale_fac2 = 1.0_prt;
                    if (up2_mag > 0.0_prt) {
                        const auto Etot2 = (Ek2 + Ek_product2) * PhysConst::q_e;
                        scale_fac2 = static_cast<amrex::ParticleReal>(
                            std::sqrt(Etot2 * (Etot2 + 2.0_prt*m2*c2) / c2) / m2
                        ) / up2_mag;
                    }
                    soa_2.m_rdata[PIdx::ux][p_pair_indices_2[i]] *= scale_fac2;
                    soa_2.m_rdata[PIdx::uy][p_pair_indices_2[i]] *= scale_fac2;
                    soa_2.m_rdata[PIdx::uz][p_pair_indices_2[i]] *= scale_fac2;
                }
            }
        });

        // Initialize the user runtime components
        for (int i = 0; i < m_num_product_species; i++)
        {
            const int start_index = int(products_np[i]);
            const int stop_index  = int(products_np[i] + num_added_vec[i]);
            ParticleCreation::DefaultInitializeRuntimeAttributes(*tile_products[i],
                                       0, 0,
                                       pc_products[i]->getUserRealAttribs(), pc_products[i]->getUserIntAttribs(),
                                       pc_products[i]->getParticleComps(), pc_products[i]->getParticleiComps(),
                                       pc_products[i]->getUserRealAttribParser(),
                                       pc_products[i]->getUserIntAttribParser(),
#ifdef WARPX_QED
                                       false, // do not initialize QED quantities, since they were initialized
                                              // when calling the SmartCopy functors
                                       pc_products[i]->get_breit_wheeler_engine_ptr(),
                                       pc_products[i]->get_quantum_sync_engine_ptr(),
#endif
                                       pc_products[i]->getIonizationInitialLevel(),
                                       start_index, stop_index);
        }

        amrex::Gpu::synchronize();
        return num_added_vec;
    }

private:
    // How many different type of species the collision produces
    int m_num_product_species;
    // Minimum reaction weight for which particles are split into new particles
    amrex::ParticleReal m_min_splitting_weight = 0.0;
    // Vectors of size m_num_product_species storing how many particles of a given species are
    // produced by a collision event. These vectors are duplicated (one version for host and one
    // for device) which is necessary with GPUs but redundant on CPU.
    amrex::Gpu::DeviceVector<int> m_num_products_device;
    amrex::Gpu::HostVector<int> m_num_products_host;
    CollisionType m_collision_type;
};
#endif // WARPX_SPLIT_AND_SCATTER_FUNC_H_
