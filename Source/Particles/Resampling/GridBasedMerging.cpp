/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "GridBasedMerging.H"

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

GridBasedMerging::GridBasedMerging (const std::string species_name)
{
    using namespace amrex::literals;

    const amrex::ParmParse pp_species_name(species_name);

    utils::parser::queryWithParser(
        pp_species_name, "resampling_algorithm_min_ppc", m_min_ppc
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

void GridBasedMerging::operator() (WarpXParIter& pti, const int lev,
                                   WarpXParticleContainer * const pc) const
{
    using namespace amrex::literals;

    auto& ptile = pc->ParticlesAt(lev, pti);
    auto& soa = ptile.GetStructOfArrays();
    auto * const AMREX_RESTRICT x = soa.GetRealData(PIdx::x).data();
    auto * const AMREX_RESTRICT y = soa.GetRealData(PIdx::y).data();
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
    auto sorted_indices_data = sorted_indices.data();

    int Ntheta = m_ntheta;
    int Nphi = m_nphi;

    auto dr = m_delta_ur;
    auto dtheta = 2.0_prt * MathConst::pi / Ntheta;
    auto dphi = MathConst::pi / Nphi;

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

    //amrex::AllPrint() << "Starting merging with " << cell_numparts << " particles in the spatial cell\n";
    //amrex::Print() << "dtheta = " << dtheta << " dphi = " << dphi << " dr = " << dr << " Ntheta = " << Ntheta << " Nphi = " << Nphi << std::endl;

            // Loop over particles and label them with the appropriate
            // momentum bin number
            for (int i = cell_start; i < cell_stop; ++i)
            {
                // get polar components of the velocity vector
                auto u_mag = std::sqrt(
                    ux[indices[i]]*ux[indices[i]]
                    + uy[indices[i]]*uy[indices[i]]
                    + uz[indices[i]]*uz[indices[i]]
                );
                auto u_theta = std::atan2(uy[indices[i]], ux[indices[i]]) + MathConst::pi;
                auto u_phi = std::acos(uz[indices[i]]/u_mag);

                int ii = static_cast<int>(u_theta / dtheta);
                int jj = static_cast<int>(u_phi / dphi);
                int kk = static_cast<int>(u_mag / dr);

                // note that the momentum bin number indexing is based on the
                // cell sorted indexing, not the particle indexing
                momentum_bin_number_data[i] = (
                    ii + jj * Ntheta + kk * Ntheta * Nphi
                );

                //amrex::Print() << "Labeling part #" << i << " : " << u_mag << ", " << u_theta << ", " << u_phi << " => " << ii << " " << jj << " " << kk << " => label = " << momentum_bin_number_data[i] << std::endl;
            }

            //amrex::AllPrint() << "done with labelling\n";

            // assign sorted_indices initial values
            // std::iota(sorted_indices.begin() + cell_start, sorted_indices.begin() + cell_stop, 0);
            for (int i = cell_start; i < cell_stop; i++)
            {
                sorted_indices_data[i] = i;
            }
            // sort indexes based on comparing values in momentum_bin_number
            // std::stable_sort(
            //     sorted_indices_data.begin() + cell_start, sorted_indices_data.begin() + cell_stop,
            //     [&momentum_bin_number_data](size_t i1, size_t i2) {
            //         return momentum_bin_number_data[i1] < momentum_bin_number_data[i2];
            //     }
            // );
            heapSort(sorted_indices_data, momentum_bin_number_data, cell_start, cell_stop);

            //amrex::AllPrint() << "done with argsort\n";

            // start by setting the running tallies equal to the first particle's attributes
            amrex::ParticleReal total_weight = w[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_x = total_weight*x[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_y = total_weight*y[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_z = total_weight*z[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_ux = total_weight*ux[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_uy = total_weight*uy[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal cluster_uz = total_weight*uz[indices[sorted_indices_data[cell_start]]];
            amrex::ParticleReal total_energy = 0.5_prt * mass * total_weight * (
                cluster_ux*cluster_ux + cluster_uy*cluster_uy + cluster_uz*cluster_uz
            );

            int particles_in_bin = 1;

            // Finally, loop through the particles in the cell and merge
            // ones in the same momentum bin
            for (int i = cell_start+1; i < cell_stop; ++i)
            {
                //amrex::Print() << "Particles in momentum cell: " << particles_in_bin << " Current particle in cell #" << momentum_bin_number_data[sorted_indices_data[i]] << std::endl;
                //amrex::Print() << "   Particle index: " << sorted_indices_data[i] << " with cell starting at " << cell_start << " and ending at " << cell_stop << std::endl;

                // check if this particle is in a new momentum bin
                if (momentum_bin_number_data[sorted_indices_data[i]] != momentum_bin_number_data[sorted_indices_data[i - 1]])
                {
                    // amrex::Print() << "Cluster merge considered\n";
                    // check if the previous bin had more than 2 particles in it
                    if (particles_in_bin > 2){
                        // get average quantities for the previous bin
                        cluster_x /= total_weight;
                        cluster_y /= total_weight;
                        cluster_z /= total_weight;
                        cluster_ux /= total_weight;
                        cluster_uy /= total_weight;
                        cluster_uz /= total_weight;

                        // perform merging of previous momentum bin particles
                        auto u_perp2 = cluster_ux*cluster_x + cluster_uy*cluster_y;
                        auto u_perp = std::sqrt(u_perp2);
                        auto cluster_u_mag2 = u_perp2 + cluster_z*cluster_z;
                        auto cluster_u_mag = std::sqrt(cluster_u_mag2);

                        // calculate required velocity magnitude to achieve
                        // energy conservation
                        auto v_mag2 = 2.0 * total_energy / (total_weight * mass);
                        auto v_perp = std::sqrt(v_mag2 - cluster_u_mag2);

                        // choose random angle for new velocity vector
                        auto phi = amrex::Random(engine) * 2.0 * MathConst::pi;

                        // set new velocity components based on chosen phi
                        auto vx = v_perp * std::cos(phi);
                        auto vy = v_perp * std::sin(phi);

                        // calculate rotation angles to parallel coord. frame
                        auto cos_theta = cluster_uz / cluster_u_mag;
                        auto sin_theta = u_perp / cluster_u_mag;
                        auto cos_phi = cluster_ux / u_perp;
                        auto sin_phi = cluster_uy / u_perp;

                        // rotate new velocity vector back to labframe
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
                        auto part_idx1 = indices[sorted_indices_data[i - 1]];
                        auto part_idx2 = indices[sorted_indices_data[i - 2]];

                        w[part_idx1] = total_weight / 2.0;
                        w[part_idx2] = total_weight / 2.0;
                        x[part_idx1] = cluster_x;
                        y[part_idx1] = cluster_y;
                        z[part_idx1] = cluster_z;
                        x[part_idx2] = cluster_x;
                        y[part_idx2] = cluster_y;
                        z[part_idx2] = cluster_z;

                        ux[part_idx1] = ux_new;
                        uz[part_idx1] = uz_new;
                        uy[part_idx1] = uy_new;
                        ux[part_idx2] = 2.0 * cluster_ux - ux_new;
                        uy[part_idx2] = 2.0 * cluster_uy - uz_new;
                        uz[part_idx2] = 2.0 * cluster_uz - uy_new;

                        // set ids of merged particles so they will be removed
                        for (int j = 2; j < particles_in_bin; j++){
                            idcpu[indices[sorted_indices_data[i - 1 - j]]] = amrex::ParticleIdCpus::Invalid;
                        }
                    }
                    // restart the tallies
                    particles_in_bin = 0;
                    total_weight = 0._prt;
                    total_energy = 0._prt;
                    cluster_x = 0_prt;
                    cluster_y = 0_prt;
                    cluster_z = 0_prt;
                    cluster_ux = 0_prt;
                    cluster_uy = 0_prt;
                    cluster_uz = 0_prt;
                }

                particles_in_bin += 1;

                cluster_x += w[indices[sorted_indices_data[i]]]*x[indices[sorted_indices_data[i]]];
                cluster_y += w[indices[sorted_indices_data[i]]]*y[indices[sorted_indices_data[i]]];
                cluster_z += w[indices[sorted_indices_data[i]]]*z[indices[sorted_indices_data[i]]];
                cluster_ux += w[indices[sorted_indices_data[i]]]*ux[indices[sorted_indices_data[i]]];
                cluster_uy += w[indices[sorted_indices_data[i]]]*uy[indices[sorted_indices_data[i]]];
                cluster_uz += w[indices[sorted_indices_data[i]]]*uz[indices[sorted_indices_data[i]]];
                total_weight += w[indices[sorted_indices_data[i]]];
                total_energy += 0.5_prt * mass * w[indices[sorted_indices_data[i]]] * (
                    ux[indices[sorted_indices_data[i]]]*ux[indices[sorted_indices_data[i]]]
                    + uy[indices[sorted_indices_data[i]]]*uy[indices[sorted_indices_data[i]]]
                    + uz[indices[sorted_indices_data[i]]]*uz[indices[sorted_indices_data[i]]]
                );
            }
        }
    );

    // amrex::Print() << "Merging complete." << std::endl;
}
