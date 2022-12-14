/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef SPLIT_AND_SCATTER_FUNC_H_
#define SPLIT_AND_SCATTER_FUNC_H_

#include "Particles/Collision/ScatteringProcess.H"

/**
 * \brief ...
 *
 * \return num_added the number of particles that were written to dst.
 */
template <typename index_type, typename Tile, typename Index,
          typename CopyFunc, typename ParticleType,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index splitScatteringParticles (
    const index_type& n_total_pairs,
    Tile& dst1, Tile& dst2, Tile& src1, Tile& src2,
    Index* mask,
    CopyFunc&& copy1, CopyFunc&& copy2,
    ParticleType* particle_ptr_1, ParticleType*  particle_ptr_2,
    const amrex::ParticleReal m1, const amrex::ParticleReal m2,
    const index_type* AMREX_RESTRICT p_pair_indices_1,
    const index_type* AMREX_RESTRICT p_pair_indices_2,
    const amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight ) noexcept
{
    using namespace amrex;

    if (n_total_pairs == 0) return 0;

    Gpu::DeviceVector<Index> offsets(n_total_pairs);
    auto p_offsets = offsets.dataPtr();
    Index offset = 0;
    for (Index i=0; i<n_total_pairs; i++)
    {
        p_offsets[i] = offset;
        if (mask[i])
        {
            offset++;
        }
        // Print() << i << "/" << n_total_pairs << " : " << offset << " " << p_offsets[i] << std::endl;
    }

    const Index num_added = mask[n_total_pairs-1] + offset;
    const auto dst1_index = dst1.numParticles();
    const auto dst2_index = dst2.numParticles();
    dst1.resize(dst1_index + num_added);
    dst2.resize(dst2_index + num_added);

    const auto src1_data = src1.getParticleTileData();
    const auto src2_data = src2.getParticleTileData();
    const auto dst1_data = dst1.getParticleTileData();
    const auto dst2_data = dst2.getParticleTileData();

    amrex::ParallelForRNG(n_total_pairs,
    [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
    {
        if (mask[i])
        {
            // Print() << p_offsets[i] << std::endl;
            // First, make copies of the colliding particles
            copy1(dst1_data, src1_data, p_pair_indices_1[i], p_offsets[i] + dst1_index, engine);
            copy2(dst2_data, src2_data, p_pair_indices_2[i], p_offsets[i] + dst2_index, engine);

            // Now we adjust the properties of the original and child particles,
            // starting with the parent particles
            auto& w1 = src1_data.m_rdata[PIdx::w][p_pair_indices_1[i]];
            auto& w2 = src2_data.m_rdata[PIdx::w][p_pair_indices_2[i]];

            // amrex::Print() << " w1 = " << w1 << std::endl;

            // Remove p_pair_reaction_weight[i] from the colliding particles' weights
            w1 -= p_pair_reaction_weight[i];
            w2 -= p_pair_reaction_weight[i];

            // If the colliding particle weight decreases to zero, remove particle by
            // setting its id to -1
            constexpr amrex::Long minus_one_long = -1;
            if (w1 <= 0.001_prt)
            {
                particle_ptr_1[p_pair_indices_1[i]].atomicSetID(minus_one_long);
            }
            if (w2 <= 0.001_prt)
            {
                particle_ptr_2[p_pair_indices_2[i]].atomicSetID(minus_one_long);
            }

            // Set the child particle properties appropriately
            dst1_data.m_rdata[PIdx::w][p_offsets[i] + dst1_index] = p_pair_reaction_weight[i];
            dst2_data.m_rdata[PIdx::w][p_offsets[i] + dst2_index] = p_pair_reaction_weight[i];

            auto& ux1 = dst1_data.m_rdata[PIdx::ux][p_offsets[i] + dst1_index];
            auto& uy1 = dst1_data.m_rdata[PIdx::uy][p_offsets[i] + dst1_index];
            auto& uz1 = dst1_data.m_rdata[PIdx::uz][p_offsets[i] + dst1_index];
            auto& ux2 = dst2_data.m_rdata[PIdx::ux][p_offsets[i] + dst2_index];
            auto& uy2 = dst2_data.m_rdata[PIdx::uy][p_offsets[i] + dst2_index];
            auto& uz2 = dst2_data.m_rdata[PIdx::uz][p_offsets[i] + dst2_index];

            // amrex::Print() << " w1 = " << w1 << "  w_new " << dst1_data.m_rdata[PIdx::w][p_offsets[i] + dst1_index] << std::endl;
            // amrex::Print() << "    -> " << std::sqrt(ux2*ux2+uy2*uy2+uz2*uz2) << "  " << src2_data.m_rdata[PIdx::w][p_pair_indices_2[i]] << "\n\n";
            // auto const vel_std = sqrt(PhysConst::kb * 300.0 / 6.67e-27);
            // auto ux2 = vel_std * amrex::RandomNormal(0_prt, 1.0_prt, engine);
            // auto uy2 = vel_std * amrex::RandomNormal(0_prt, 1.0_prt, engine);
            // auto uz2 = vel_std * amrex::RandomNormal(0_prt, 1.0_prt, engine);


            // for simplicity and for now we assume non-relativistic particles
            auto const uCOM_x = (m1 * ux1 + m2 * ux2) / (m1 + m2);
            auto const uCOM_y = (m1 * uy1 + m2 * uy2) / (m1 + m2);
            auto const uCOM_z = (m1 * uz1 + m2 * uz2) / (m1 + m2);

            // transform to COM frame
            ux1 -= uCOM_x;
            uy1 -= uCOM_y;
            uz1 -= uCOM_z;
            // ux2 -= uCOM_x;
            // uy2 -= uCOM_y;
            // uz2 -= uCOM_z;

            if (mask[i] == int(ScatteringProcessType::ELASTIC)) {
                // randomly rotate the velocity vector for the first particle
                ParticleUtils::RandomizeVelocity(
                    ux1, uy1, uz1, sqrt(ux1*ux1 + uy1*uy1 + uz1*uz1), engine
                );
                // set the second particles velocity so that the total momentum
                // is zero
                // ux2 = -ux1 * m1 / m2;
                // uy2 = -uy1 * m1 / m2;
                // uz2 = -uz1 * m1 / m2;
            } else if (mask[i] == int(ScatteringProcessType::BACK)) {
                // reverse the velocity vectors of both particles
                ux1 *= -1.0_prt;
                uy1 *= -1.0_prt;
                uz1 *= -1.0_prt;
                // ux2 *= -1.0_prt;
                // uy2 *= -1.0_prt;
                // uz2 *= -1.0_prt;
            } else if (mask[i] == int(ScatteringProcessType::CHARGE_EXCHANGE)) {
                if (m1 == m2) {
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
                    Abort("Uneven mass charge-exchange not implemented yet.");
                }
            }
            else {
                Abort("Unknown scattering process.");
            }
            // transform back to labframe
            ux1 += uCOM_x;
            uy1 += uCOM_y;
            uz1 += uCOM_z;
            // ux2 += uCOM_x;
            // uy2 += uCOM_y;
            // uz2 += uCOM_z;
        }
    });

    Gpu::synchronize();
    return num_added;
}

#endif // SPLIT_AND_SCATTER_FUNC_H_