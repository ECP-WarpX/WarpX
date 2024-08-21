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
        "Collider-relevant diagnostics must involve exactly two species");

    // RZ coordinate is not supported
#if (defined WARPX_DIM_RZ)
    WARPX_ABORT_WITH_MESSAGE(
        "Collider-relevant diagnostics do not work in RZ geometry.");
#endif

    ablastr::warn_manager::WMRecordWarning(
        "Diagnostics",
        "The collider-relevant reduced diagnostic is meant for \
        colliding species propagating along the z direction.",
        ablastr::warn_manager::WarnPriority::low);

    ablastr::warn_manager::WMRecordWarning(
        "Diagnostics",
        "The collider-relevant reduced diagnostic only considers the \
        coarsest level of refinement for the calculations involving chi.",
        ablastr::warn_manager::WarnPriority::low);

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

    // get WarpX class object
    auto& warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    const MultiParticleContainer& mypc =  warpx.GetPartContainer();

    // resize data array
    m_data.resize(m_bin_num,0.0_rt);


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
                ofs << "bin" + std::to_string(1+i)
                             + "=" + std::to_string(b) + "()";
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

    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get a reference to WarpX instance
    auto& warpx = WarpX::GetInstance();

    // get time at level 0
    auto const t = warpx.gett_new(0);

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

    // zero-out old data on the host
    std::fill(m_data.begin(), m_data.end(), Real(0.0));
    Gpu::DeviceVector< Real > d_data( m_data.size(), 0.0 );
    Real* const AMREX_RESTRICT dptr_data = d_data.dataPtr();

    // Enable tiling
    amrex::MFItInfo info;
    if (amrex::Gpu::notInLaunchRegion()) { info.EnableTiling(species_1.tile_size); }

    int const nlevs = std::max(0, species_1.finestLevel()+1); // species_1 ?
    for (int lev = 0; lev < nlevs; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

        constexpr int getpos_offset = 0; // NOT SURE
        for (amrex::MFIter mfi = species_1.MakeMFIter(lev, info); mfi.isValid(); ++mfi){

            ParticleTileType& ptile_1 = species_1.ParticlesAt(lev, mfi);
            ParticleTileType& ptile_2 = species_2.ParticlesAt(lev, mfi);

            ParticleBins bins_1 = ParticleUtils::findParticlesInEachCell( lev, mfi, ptile_1 );
            ParticleBins bins_2 = ParticleUtils::findParticlesInEachCell( lev, mfi, ptile_2 );

            // Extract low-level data
            auto const n_cells = static_cast<int>(bins_1.numBins());

            // - Species 1
            const auto soa_1 = ptile_1.getParticleTileData();
            const auto soa_2 = ptile_2.getParticleTileData();


            // Loop over cells
            amrex::ParallelFor( n_cells,
                [=] AMREX_GPU_DEVICE (int i_cell) noexcept
            {
                amrex::ParticleReal * const AMREX_RESTRICT w1  = soa_1.m_rdata[PIdx::w];
                amrex::ParticleReal * const AMREX_RESTRICT u1x = soa_1.m_rdata[PIdx::ux];
                amrex::ParticleReal * const AMREX_RESTRICT u1y = soa_1.m_rdata[PIdx::uy];
                amrex::ParticleReal * const AMREX_RESTRICT u1z = soa_1.m_rdata[PIdx::uz];

                amrex::ParticleReal * const AMREX_RESTRICT w2  = soa_2.m_rdata[PIdx::w];
                amrex::ParticleReal * const AMREX_RESTRICT u2x = soa_2.m_rdata[PIdx::ux];
                amrex::ParticleReal * const AMREX_RESTRICT u2y = soa_2.m_rdata[PIdx::uy];
                amrex::ParticleReal * const AMREX_RESTRICT u2z = soa_2.m_rdata[PIdx::uz];




                amrex::Print(0) << soa_1.m_size << " +++++++++ " << soa_2.m_size << std::endl;

                for(int i_1=0; i_1<soa_1.m_size; i_1++){
                    for(int i_2=0; i_2<soa_2.m_size; i_2++){

                        Real m1_square = m1*m1*w1[i_1]*w1[i_1];
                        Real m2_square = m2*m2*w2[i_2]*w2[i_2];
                        Real gamma1 = std::sqrt(1. + u1x[i_1]*u1x[i_1] + u1y[i_1]*u1y[i_1] + u1z[i_1]*u1z[i_1]);
                        Real gamma2 = std::sqrt(1. + u2x[i_2]*u2x[i_2] + u2y[i_2]*u2y[i_2] + u2z[i_2]*u2z[i_2]);
                        Real u1_dot_u2 = u1x[i_1]*u2x[i_2] + u1y[i_1]*u2y[i_2] + u1z[i_1]*u2z[i_2];
                        //Real E_com = PhysConst::c*PhysConst::c *std::sqrt(m1_square + m2_square + 2*m1*w1[i_1]*m2*w2[i_2] * (gamma1*gamma2 - u1_dot_u2)); // cancellation errors?
                        Real E_com = PhysConst::c*PhysConst::c *std::sqrt(m1_square + m2_square + 2*m1*w1[i_1]*m2*w2[i_2] * (gamma1*gamma2 - u1_dot_u2)); // cancellation errors?
                        // determine particle bin
                        int const bin = int(Math::floor((E_com-bin_min)/bin_size));
                        if ( bin<0 || bin>=num_bins ) { return; } // discard if out-of-range

                        // add value to dL_dEcom_dt bin
                        amrex::HostDevice::Atomic::Add(&dptr_data[bin], w1[i_1]*w2[i_2]);

                        amrex::Print(0) << w1[i_1] << " ooo " << w2[i_2] << std::endl;
                    } // particles species 2
                } // particles species 1
            }); // cells
        } // boxes
    } // levels

    // blocking copy from device to host
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
        d_data.begin(), d_data.end(), m_data.begin());

    // reduced sum over mpi ranks
    ParallelDescriptor::ReduceRealSum
        (m_data.data(), static_cast<int>(m_data.size()), ParallelDescriptor::IOProcessorNumber());

    // normalize data
    for ( int i = 0; i < m_bin_num; ++i )
    {
        m_data[i] *= 2._rt * PhysConst::c;
    }

#endif // not RZ
}
