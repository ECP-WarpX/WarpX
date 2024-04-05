/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "VelocityCoincidenceThinning.H"

#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/ParticleUtils.H"
#include "Utils/TextMsg.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_DenseBins.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_Particles.H>
#include <AMReX_Random.H>
#include <AMReX_StructOfArrays.H>

#include <AMReX_BaseFwd.H>

VelocityCoincidenceThinning::VelocityCoincidenceThinning (const std::string& species_name)
{
    using namespace amrex::literals;

    const amrex::ParmParse pp_species_name(species_name);

    utils::parser::queryWithParser(
        pp_species_name, "resampling_min_ppc", m_min_ppc
    );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_min_ppc >= 1,
        "Resampling min_ppc should be greater than or equal to 1"
    );

    utils::parser::queryWithParser(
        pp_species_name, "resampling_algorithm_delta_ur", m_delta_ur
    );
    utils::parser::queryWithParser(
        pp_species_name, "resampling_algorithm_n_theta", m_ntheta
    );
    utils::parser::queryWithParser(
        pp_species_name, "resampling_algorithm_n_phi", m_nphi
    );
}

void VelocityCoincidenceThinning::operator() (WarpXParIter& pti, const int lev,
                                   WarpXParticleContainer * const pc) const
{
    using namespace amrex::literals;

    auto& ptile = pc->ParticlesAt(lev, pti);
    auto& soa = ptile.GetStructOfArrays();
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D)
    auto * const AMREX_RESTRICT x = soa.GetRealData(PIdx::x).data();
#elif defined(WARPX_DIM_RZ)
    auto * const AMREX_RESTRICT x = soa.GetRealData(PIdx::x).data(); // rename to PIdx::r after PR #4667
#endif
#if defined(WARPX_DIM_3D)
    auto * const AMREX_RESTRICT y = soa.GetRealData(PIdx::y).data();
#endif
    auto * const AMREX_RESTRICT z = soa.GetRealData(PIdx::z).data();
    auto * const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
    auto * const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
    auto * const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();
    auto * const AMREX_RESTRICT w = soa.GetRealData(PIdx::w).data();
    auto * const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();

    // Using this function means that we must loop over the cells in the ParallelFor.
    auto bins = ParticleUtils::findParticlesInEachCell(lev, pti, ptile);

    const auto n_cells = static_cast<int>(bins.numBins());
    auto *const indices = bins.permutationPtr();
    auto *const cell_offsets = bins.offsetsPtr();

    const int min_ppc = m_min_ppc;

    const auto mass = pc->getMass();

    // create a GPU vector to hold the momentum cluster index for each particle
    amrex::Gpu::DeviceVector<int> momentum_bin_number;
    momentum_bin_number.resize(bins.numItems());
    auto* momentum_bin_number_data = momentum_bin_number.data();

    // create a GPU vector to hold the index sorting for the momentum bins
    amrex::Gpu::DeviceVector<int> sorted_indices;
    sorted_indices.resize(bins.numItems());
    auto* sorted_indices_data = sorted_indices.data();

    const int Ntheta = m_ntheta;
    const int Nphi = m_nphi;

    auto dr = m_delta_ur;
    auto dtheta = 2.0_prt * MathConst::pi / Ntheta;
    auto dphi = MathConst::pi / Nphi;

    auto heapSort = HeapSort();

    // Loop over cells
    amrex::ParallelForRNG( n_cells,
        [=] AMREX_GPU_DEVICE (int i_cell, amrex::RandomEngine const& engine) noexcept
        {
            // The particles that are in the cell `i_cell` are
            // given by the `indices[cell_start:cell_stop]`
            const auto cell_start = static_cast<int>(cell_offsets[i_cell]);
            const auto cell_stop  = static_cast<int>(cell_offsets[i_cell+1]);
            const int cell_numparts = cell_stop - cell_start;

            // do nothing for cells with less particles than min_ppc
            // (this intentionally includes skipping empty cells, too)
            if (cell_numparts < min_ppc) {
                return;
            }

            // Loop over particles and label them with the appropriate
            // momentum bin number
            for (int i = cell_start; i < cell_stop; ++i)
            {
                // get polar components of the velocity vector
                auto u_mag = std::sqrt(
                    ux[indices[i]]*ux[indices[i]] +
                    uy[indices[i]]*uy[indices[i]] +
                    uz[indices[i]]*uz[indices[i]]
                );
                auto u_theta = std::atan2(uy[indices[i]], ux[indices[i]]) + MathConst::pi;
                auto u_phi = std::acos(uz[indices[i]]/u_mag);

                const int ii = static_cast<int>(u_theta / dtheta);
                const int jj = static_cast<int>(u_phi / dphi);
                const int kk = static_cast<int>(u_mag / dr);

                // note that the momentum bin number indexing is based on the
                // cell sorted indexing, not the particle indexing
                momentum_bin_number_data[i] = (
                    ii + jj * Ntheta + kk * Ntheta * Nphi
                );
            }

            // assign sorted_indices initial values
            for (int i = cell_start; i < cell_stop; i++)
            {
                sorted_indices_data[i] = i;
            }
            // sort indexes based on comparing values in momentum_bin_number
            heapSort(sorted_indices_data, momentum_bin_number_data, cell_start, cell_numparts);

            // start by setting the running tallies equal to the first particle's attributes
            amrex::ParticleReal total_weight = w[indices[sorted_indices_data[cell_start]]];
#if !defined(WARPX_DIM_1D_Z)
            amrex::ParticleReal cluster_x = total_weight*x[indices[sorted_indices_data[cell_start]]];
#endif
#if defined(WARPX_DIM_3D)
            amrex::ParticleReal cluster_y = total_weight*y[indices[sorted_indices_data[cell_start]]];
#endif
            amrex::ParticleReal cluster_z = total_weight*z[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_ux = total_weight*ux[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_uy = total_weight*uy[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_uz = total_weight*uz[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal total_energy = 0.5_prt * mass * total_weight * (
                ux[indices[sorted_indices_data[cell_start]]]*ux[indices[sorted_indices_data[cell_start]]] +
                uy[indices[sorted_indices_data[cell_start]]]*uy[indices[sorted_indices_data[cell_start]]] +
                uz[indices[sorted_indices_data[cell_start]]]*uz[indices[sorted_indices_data[cell_start]]]
            );

            int particles_in_bin = 1;

            // Finally, loop through the particles in the cell and merge
            // ones in the same momentum bin
            for (int i = cell_start+1; i < cell_stop; ++i)
            {
                // check if this particle is in a new momentum bin
                if (momentum_bin_number_data[sorted_indices_data[i]] != momentum_bin_number_data[sorted_indices_data[i - 1]])
                {
                    // check if the previous bin had more than 2 particles in it
                    if (particles_in_bin > 2){
                        // get average quantities for the previous bin
#if !defined(WARPX_DIM_1D_Z)
                        cluster_x /= total_weight;
#endif
#if defined(WARPX_DIM_3D)
                        cluster_y /= total_weight;
#endif
                        cluster_z /= total_weight;
                        cluster_ux /= total_weight;
                        cluster_uy /= total_weight;
                        cluster_uz /= total_weight;

                        // perform merging of momentum bin particles
                        auto u_perp2 = cluster_ux*cluster_ux + cluster_uy*cluster_uy;
                        auto u_perp = std::sqrt(u_perp2);
                        auto cluster_u_mag2 = u_perp2 + cluster_uz*cluster_uz;
                        auto cluster_u_mag = std::sqrt(cluster_u_mag2);

                        // calculate required velocity magnitude to achieve
                        // energy conservation
                        auto v_mag2 = 2._prt * total_energy / (total_weight * mass);
                        auto v_perp = std::sqrt(v_mag2 - cluster_u_mag2);

                        // choose random angle for new velocity vector
                        auto phi = amrex::Random(engine) * MathConst::pi;

                        // set new velocity components based on chosen phi
                        auto vx = v_perp * std::cos(phi);
                        auto vy = v_perp * std::sin(phi);

                        // calculate rotation angles to parallel coord. frame
                        auto cos_theta = cluster_uz / cluster_u_mag;
                        auto sin_theta = u_perp / cluster_u_mag;
                        auto cos_phi = cluster_ux / u_perp;
                        auto sin_phi = cluster_uy / u_perp;

                        // rotate new velocity vector to labframe
                        auto ux_new = (
                            vx * cos_theta * cos_phi - vy * sin_phi
                            + cluster_u_mag * sin_theta * cos_phi
                        );
                        auto uy_new = (
                            vx * cos_theta * sin_phi + vy * cos_phi
                            + cluster_u_mag * sin_theta * sin_phi
                        );
                        auto uz_new = -vx * sin_theta + cluster_u_mag * cos_theta;

                        // set the previous two particles' attributes according to
                        // the previous bin's values
                        const auto part_idx1 = indices[sorted_indices_data[i - 1]];
                        const auto part_idx2 = indices[sorted_indices_data[i - 2]];

                        w[part_idx1] = total_weight / 2._prt;
                        w[part_idx2] = total_weight / 2._prt;
#if !defined(WARPX_DIM_1D_Z)
                        x[part_idx1] = cluster_x;
                        x[part_idx2] = cluster_x;
#endif
#if defined(WARPX_DIM_3D)
                        y[part_idx1] = cluster_y;
                        y[part_idx2] = cluster_y;
#endif
                        z[part_idx1] = cluster_z;
                        z[part_idx2] = cluster_z;

                        ux[part_idx1] = ux_new;
                        uy[part_idx1] = uy_new;
                        uz[part_idx1] = uz_new;
                        ux[part_idx2] = 2._prt * cluster_ux - ux_new;
                        uy[part_idx2] = 2._prt * cluster_uy - uy_new;
                        uz[part_idx2] = 2._prt * cluster_uz - uz_new;

                        // set ids of merged particles so they will be removed
                        for (int j = 2; j < particles_in_bin; j++){
                            idcpu[indices[sorted_indices_data[i - 1 - j]]] = amrex::ParticleIdCpus::Invalid;
                        }
                    }
                    // restart the tallies
                    particles_in_bin = 0;
                    total_weight = 0._prt;
                    total_energy = 0._prt;
#if !defined(WARPX_DIM_1D_Z)
                    cluster_x = 0_prt;
#endif
#if defined(WARPX_DIM_3D)
                    cluster_y = 0_prt;
#endif
                    cluster_z = 0_prt;
                    cluster_ux = 0_prt;
                    cluster_uy = 0_prt;
                    cluster_uz = 0_prt;
                }

                particles_in_bin += 1;
                const auto part_idx = indices[sorted_indices_data[i]];

#if !defined(WARPX_DIM_1D_Z)
                cluster_x += w[part_idx]*x[part_idx];
#endif
#if defined(WARPX_DIM_3D)
                cluster_y += w[part_idx]*y[part_idx];
#endif
                cluster_z += w[part_idx]*z[part_idx];
                cluster_ux += w[part_idx]*ux[part_idx];
                cluster_uy += w[part_idx]*uy[part_idx];
                cluster_uz += w[part_idx]*uz[part_idx];
                total_weight += w[part_idx];
                total_energy += 0.5_prt * mass * w[part_idx] * (
                    ux[part_idx]*ux[part_idx] + uy[part_idx]*uy[part_idx]
                    + uz[part_idx]*uz[part_idx]
                );
            }
        }
    );
}
