/* Copyright 2023 The WarpX Community
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
#include "Particles/NamedComponentParticleContainer.H"

#include <AMReX_Particle.H>

/**
 * \brief Function that performs the particle scattering and injection due
 * to binary collisions.
 *
 * \return num_added the number of particles added to each species.
 */
template <typename index_type, typename Tile, typename Index, typename CopyFunc,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
int splitScatteringParticles (
    const index_type& n_total_pairs,
    Tile& ptile1, Tile& ptile2,
    const Index* AMREX_RESTRICT mask,
    CopyFunc&& copy1, CopyFunc&& copy2,
    const amrex::ParticleReal m1, const amrex::ParticleReal m2,
    const index_type* AMREX_RESTRICT p_pair_indices_1,
    const index_type* AMREX_RESTRICT p_pair_indices_2,
    const amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight ) noexcept
{
    using namespace amrex;

    if (n_total_pairs == 0) {
        return 0;
    }

    Gpu::DeviceVector<Index> offsets(n_total_pairs);
    Index* AMREX_RESTRICT p_offsets = offsets.dataPtr();

    // The following is used to calculate the appropriate offsets. Note that
    // a standard cummulative sum is not appropriate since the mask is also
    // used to specify the type of collision and can therefore have values >1
    auto const num_added = amrex::Scan::PrefixSum<Index>(n_total_pairs,
        [=] AMREX_GPU_DEVICE (Index i) -> Index { return mask[i] ? 1 : 0; },
        [=] AMREX_GPU_DEVICE (Index i, Index s) { p_offsets[i] = s; },
        amrex::Scan::Type::exclusive, amrex::Scan::retSum
    );

    const auto ptile1_index = ptile1.numParticles();
    const auto ptile2_index = ptile2.numParticles();
    ptile1.resize(ptile1_index + num_added);
    ptile2.resize(ptile2_index + num_added);

    const auto ptile1_data = ptile1.getParticleTileData();
    const auto ptile2_data = ptile2.getParticleTileData();

    ParallelForRNG(n_total_pairs,
    [=] AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
    {
        if (mask[i])
        {
            // First, make copies of the colliding particles
            copy1(ptile1_data, ptile1_data, p_pair_indices_1[i], p_offsets[i] + ptile1_index, engine);
            copy2(ptile2_data, ptile2_data, p_pair_indices_2[i], p_offsets[i] + ptile2_index, engine);

            // Now we adjust the properties of the original and child particles,
            // starting with the parent particles
            auto& w1 = ptile1_data.m_rdata[PIdx::w][p_pair_indices_1[i]];
            auto& w2 = ptile2_data.m_rdata[PIdx::w][p_pair_indices_2[i]];
            uint64_t* AMREX_RESTRICT idcpu1 = ptile1_data.m_idcpu;
            uint64_t* AMREX_RESTRICT idcpu2 = ptile2_data.m_idcpu;

            // Note: Particle::atomicSetID should also be provided as a standalone helper function in AMReX
            //       to replace the following lambda.
            auto const atomicSetIdMinus = [] AMREX_GPU_DEVICE (uint64_t & idcpu)
            {
#if defined(AMREX_USE_OMP)
#pragma omp atomic write
                idcpu = amrex::ParticleIdCpus::Invalid;
#else
                amrex::Gpu::Atomic::Exch(
                    (unsigned long long *)&idcpu,
                    (unsigned long long)amrex::ParticleIdCpus::Invalid
                );
#endif
            };

            // Remove p_pair_reaction_weight[i] from the colliding particles' weights.
            // If the colliding particle weight decreases to zero, remove particle by
            // setting its id to -1.
            Gpu::Atomic::AddNoRet(&w1, -p_pair_reaction_weight[i]);
            if (w1 <= 0._prt) {
                atomicSetIdMinus(idcpu1[p_pair_indices_1[i]]);
            }

            Gpu::Atomic::AddNoRet(&w2, -p_pair_reaction_weight[i]);
            if (w2 <= 0._prt) {
                atomicSetIdMinus(idcpu2[p_pair_indices_2[i]]);
            }

            // Set the child particle properties appropriately
            ptile1_data.m_rdata[PIdx::w][p_offsets[i] + ptile1_index] = p_pair_reaction_weight[i];
            ptile2_data.m_rdata[PIdx::w][p_offsets[i] + ptile2_index] = p_pair_reaction_weight[i];

            auto& ux1 = ptile1_data.m_rdata[PIdx::ux][p_offsets[i] + ptile1_index];
            auto& uy1 = ptile1_data.m_rdata[PIdx::uy][p_offsets[i] + ptile1_index];
            auto& uz1 = ptile1_data.m_rdata[PIdx::uz][p_offsets[i] + ptile1_index];
            auto& ux2 = ptile2_data.m_rdata[PIdx::ux][p_offsets[i] + ptile2_index];
            auto& uy2 = ptile2_data.m_rdata[PIdx::uy][p_offsets[i] + ptile2_index];
            auto& uz2 = ptile2_data.m_rdata[PIdx::uz][p_offsets[i] + ptile2_index];

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
            ux2 += uCOM_x;
            uy2 += uCOM_y;
            uz2 += uCOM_z;
        }
    });

    Gpu::synchronize();
    return static_cast<int>(num_added);
}
#endif // SPLIT_AND_SCATTER_FUNC_H_
