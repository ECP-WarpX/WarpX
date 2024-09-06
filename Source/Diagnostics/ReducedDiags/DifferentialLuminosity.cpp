/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Arianna Formenti, Yinjian Zhao
 * License: BSD-3-Clause-LBNL
 */
#include "DifferentialLuminosity.H"

#include "Diagnostics/ReducedDiags/ReducedDiags.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/SpeciesPhysicalProperties.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/ParticleUtils.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Dim3.H>
#include <AMReX_Extension.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>
#include <AMReX_Vector.H>

#include <ablastr/warn_manager/WarnManager.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

using ParticleType = WarpXParticleContainer::ParticleType;
using ParticleTileType = WarpXParticleContainer::ParticleTileType;
using ParticleTileDataType = ParticleTileType::ParticleTileDataType;
using ParticleBins = amrex::DenseBins<ParticleTileDataType>;
using index_type = ParticleBins::index_type;

using namespace amrex;

DifferentialLuminosity::DifferentialLuminosity (const std::string& rd_name)
: ReducedDiags{rd_name}
{
    // read colliding species names - must be 2
    const amrex::ParmParse pp_rd_name(m_rd_name);
    pp_rd_name.getarr("species", m_beam_name);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_beam_name.size() == 2u,
        "DifferentialLuminosity diagnostics must involve exactly two species");

    // RZ coordinate is not supported
#if (defined WARPX_DIM_RZ)
    WARPX_ABORT_WITH_MESSAGE(
        "DifferentialLuminosity diagnostics do not work in RZ geometry.");
#endif

    // read bin parameters
    int bin_num = 0;
    amrex::Real bin_max = 0.0_rt, bin_min = 0.0_rt;
    utils::parser::getWithParser(pp_rd_name, "bin_number", bin_num);
    utils::parser::getWithParser(pp_rd_name, "bin_max",    bin_max);
    utils::parser::getWithParser(pp_rd_name, "bin_min",    bin_min);
    m_bin_num = bin_num;
    m_bin_max = bin_max;
    m_bin_min = bin_min;
    m_bin_size = (bin_max - bin_min) / bin_num;

    // resize and zero-out data array
    m_data.resize(m_bin_num,0.0_rt);
    d_data.resize(m_bin_num,0.0_rt);

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if ( m_write_header )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app);
            // write header row
            int off = 0;
            ofs << "#";
            ofs << "[" << off++ << "]step()";
            ofs << m_sep;
            ofs << "[" << off++ << "]time(s)";
            for (int i = 0; i < m_bin_num; ++i)
            {
                ofs << m_sep;
                ofs << "[" << off++ << "]";
                const Real b = m_bin_min + m_bin_size*(Real(i)+0.5_rt);
                ofs << "bin" << 1+i << "=" << b << "(eV)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }
}

void DifferentialLuminosity::ComputeDiags (int step)
{
#if defined(WARPX_DIM_RZ)
    amrex::ignore_unused(step);
#else
    WARPX_PROFILE("DifferentialLuminosity::ComputeDiags");

    // Since this diagnostic *accumulates* the luminosity in the
    // array d_data, we add contributions at *each timestep*, but
    // we only write the data to file at intervals specified by the user.
    const Real c2 = PhysConst::c*PhysConst::c;
    const Real c_over_qe = PhysConst::c/PhysConst::q_e;

    // get a reference to WarpX instance
    auto& warpx = WarpX::GetInstance();
    const Real dt = warpx.getdt(0);
    // get cell volume
    Geometry const & geom = warpx.Geom(0);
    const Real dV = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));

    // declare local variables
    auto const num_bins = m_bin_num;
    Real const bin_min  = m_bin_min;
    Real const bin_size = m_bin_size;

    // get MultiParticleContainer class object
    const MultiParticleContainer& mypc = warpx.GetPartContainer();

    auto& species_1 = mypc.GetParticleContainerFromName(m_beam_name[0]);
    auto& species_2 = mypc.GetParticleContainerFromName(m_beam_name[1]);

    const ParticleReal m1 = species_1.getMass();
    const ParticleReal m2 = species_2.getMass();

    amrex::Real* const AMREX_RESTRICT dptr_data = d_data.dataPtr();

    // Enable tiling
    amrex::MFItInfo info;
    if (amrex::Gpu::notInLaunchRegion()) { info.EnableTiling(WarpXParticleContainer::tile_size); }

    int const nlevs = std::max(0, species_1.finestLevel()+1); // species_1 ?
    for (int lev = 0; lev < nlevs; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

        for (amrex::MFIter mfi = species_1.MakeMFIter(lev, info); mfi.isValid(); ++mfi){

            ParticleTileType& ptile_1 = species_1.ParticlesAt(lev, mfi);
            ParticleTileType& ptile_2 = species_2.ParticlesAt(lev, mfi);

            ParticleBins bins_1 = ParticleUtils::findParticlesInEachCell( lev, mfi, ptile_1 );
            ParticleBins bins_2 = ParticleUtils::findParticlesInEachCell( lev, mfi, ptile_2 );

            // Species
            const auto soa_1 = ptile_1.getParticleTileData();
            index_type* AMREX_RESTRICT indices_1 = bins_1.permutationPtr();
            index_type const* AMREX_RESTRICT cell_offsets_1 = bins_1.offsetsPtr();

            // Particle data in the tile/box
            amrex::ParticleReal * const AMREX_RESTRICT w1  = soa_1.m_rdata[PIdx::w];
            amrex::ParticleReal * const AMREX_RESTRICT u1x = soa_1.m_rdata[PIdx::ux];
            amrex::ParticleReal * const AMREX_RESTRICT u1y = soa_1.m_rdata[PIdx::uy]; // v*gamma=p/m
            amrex::ParticleReal * const AMREX_RESTRICT u1z = soa_1.m_rdata[PIdx::uz];

            const auto soa_2 = ptile_2.getParticleTileData();
            index_type* AMREX_RESTRICT indices_2 = bins_2.permutationPtr();
            index_type const* AMREX_RESTRICT cell_offsets_2 = bins_2.offsetsPtr();

            amrex::ParticleReal * const AMREX_RESTRICT w2  = soa_2.m_rdata[PIdx::w];
            amrex::ParticleReal * const AMREX_RESTRICT u2x = soa_2.m_rdata[PIdx::ux];
            amrex::ParticleReal * const AMREX_RESTRICT u2y = soa_2.m_rdata[PIdx::uy];
            amrex::ParticleReal * const AMREX_RESTRICT u2z = soa_2.m_rdata[PIdx::uz];

            // Extract low-level data
            auto const n_cells = static_cast<int>(bins_1.numBins());

            // Loop over cells
            amrex::ParallelFor( n_cells,
                [=] AMREX_GPU_DEVICE (int i_cell) noexcept
            {
                // The particles from species1 that are in the cell `i_cell` are
                // given by the `indices_1[cell_start_1:cell_stop_1]`
                index_type const cell_start_1 = cell_offsets_1[i_cell];
                index_type const cell_stop_1  = cell_offsets_1[i_cell+1];
                // Same for species 2
                index_type const cell_start_2 = cell_offsets_2[i_cell];
                index_type const cell_stop_2  = cell_offsets_2[i_cell+1];

                for(index_type i_1=cell_start_1; i_1<cell_stop_1; ++i_1){
                    for(index_type i_2=cell_start_2; i_2<cell_stop_2; ++i_2){

                        index_type const j_1 = indices_1[i_1];
                        index_type const j_2 = indices_2[i_2];

                        Real p1t=0, p1x=0, p1y=0, p1z=0; // components of 4-momentum of particle 1
                        Real const u1_square =  u1x[j_1]*u1x[j_1] + u1y[j_1]*u1y[j_1] + u1z[j_1]*u1z[j_1];
                        if (m1 != 0) {
                            p1t = m1*std::sqrt( c2 + u1_square );
                            p1x = m1*u1x[j_1];
                            p1y = m1*u1y[j_1];
                            p1z = m1*u1z[j_1];
                        } else { // photon case (momentum is normalized by m_e)
                            p1t = PhysConst::m_e*std::sqrt( u1_square );
                            p1x = PhysConst::m_e*u1x[j_1];
                            p1y = PhysConst::m_e*u1y[j_1];
                            p1z = PhysConst::m_e*u1z[j_1];
                        }

                        Real p2t=0, p2x=0, p2y=0, p2z=0; // components of 4-momentum of particle 2
                        Real const u2_square =  u2x[j_2]*u2x[j_2] + u2y[j_2]*u2y[j_2] + u2z[j_2]*u2z[j_2];
                        if (m2 != 0) {
                            p2t = m2*std::sqrt( c2 + u2_square );
                            p2x = m2*u2x[j_2];
                            p2y = m2*u2y[j_2];
                            p2z = m2*u2z[j_2];
                        } else { // photon case (momentum is normalized by m_e)
                            p2t = PhysConst::m_e*std::sqrt(u2_square);
                            p2x = PhysConst::m_e*u2x[j_2];
                            p2y = PhysConst::m_e*u2y[j_2];
                            p2z = PhysConst::m_e*u2z[j_2];
                        }

                        // center of mass energy in eV
                        Real const E_com = c_over_qe * std::sqrt(m1*m1*c2 + m2*m2*c2 + 2*(p1t*p2t - p1x*p2x - p1y*p2y - p1z*p2z));

                        // determine particle bin
                        int const bin = int(Math::floor((E_com-bin_min)/bin_size));

                        if ( bin<0 || bin>=num_bins ) { continue; } // discard if out-of-range

                        Real const inv_p1t = 1.0_rt/p1t;
                        Real const inv_p2t = 1.0_rt/p2t;

                        Real const beta1_minus_beta2_x = p1x * inv_p1t - p2x * inv_p2t;
                        Real const beta1_minus_beta2_y = p1y * inv_p1t - p2y * inv_p2t;
                        Real const beta1_minus_beta2_z = p1z * inv_p1t - p2z * inv_p2t;
                        Real const beta1_minus_beta2_square = beta1_minus_beta2_x*beta1_minus_beta2_x + beta1_minus_beta2_y*beta1_minus_beta2_y + beta1_minus_beta2_z*beta1_minus_beta2_z;

                        Real const beta1_cross_beta2_x = (p1y*p2z - p1z*p2y) * inv_p1t*inv_p2t;
                        Real const beta1_cross_beta2_y = (p1z*p2x - p1x*p2z) * inv_p1t*inv_p2t;
                        Real const beta1_cross_beta2_z = (p1x*p2y - p1y*p2x) * inv_p1t*inv_p2t;

                        Real const beta1_cross_beta2_square = (beta1_cross_beta2_x*beta1_cross_beta2_x + beta1_cross_beta2_y*beta1_cross_beta2_y + beta1_cross_beta2_z*beta1_cross_beta2_z);

                        Real const radicand = beta1_minus_beta2_square - beta1_cross_beta2_square;

                        Real const dL_dEcom = PhysConst::c * std::sqrt( radicand ) * w1[j_1] * w2[j_2] / dV / bin_size * dt; // m^-2 eV^-1

                        amrex::HostDevice::Atomic::Add(&dptr_data[bin], dL_dEcom);

                    } // particles species 2
                } // particles species 1
            }); // cells
        } // boxes
    } // levels

    // Only write to file at intervals specified by the user.
    // At these intervals, the data needs to ready on the CPU,
    // so we copy it from the GPU to the CPU and reduce across MPI ranks.
    if (m_intervals.contains(step+1)) {
        // blocking copy from device to host
        amrex::Gpu::copy(amrex::Gpu::deviceToHost,
            d_data.begin(), d_data.end(), m_data.begin());

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum
            (m_data.data(), static_cast<int>(m_data.size()), ParallelDescriptor::IOProcessorNumber());
    }

#endif // not RZ

}
