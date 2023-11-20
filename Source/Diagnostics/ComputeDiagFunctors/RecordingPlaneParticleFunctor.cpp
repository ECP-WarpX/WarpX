/* Copyright 2023 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "RecordingPlaneParticleFunctor.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BaseFwd.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleContainerBase.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_ParticleTile.H>


RecordParticles::RecordParticles (const WarpXParIter &a_pti, TmpParticles& tmp_particle_data,
                                  amrex::Real z_location, int a_offset)
    : m_z_location(z_location)
{
    m_get_position = GetParticlePosition<PIdx>(a_pti, a_offset);

    const auto lev = a_pti.GetLevel();
    const auto index = a_pti.GetPairIndex();

    zp_old = tmp_particle_data[lev][index][TmpIdx::zold].dataPtr();
}

PlaneCrossingTime::PlaneCrossingTime (const WarpXParIter& a_pti, TmpParticles& tmp_particle_data,
                                      amrex::Real current_time,
                                      amrex::Real z_station_location, int a_index, int a_offset)
    : m_z_station(z_station_location), m_current_time(current_time), m_index(a_index)
{
    using namespace amrex::literals;
    m_get_position = GetParticlePosition<PIdx>(a_pti, a_offset);

    auto& attribs = a_pti.GetAttribs();
    m_wpnew = attribs[PIdx::w].dataPtr();
    m_uxnew = attribs[PIdx::ux].dataPtr();
    m_uynew = attribs[PIdx::uy].dataPtr();
    m_uznew = attribs[PIdx::uz].dataPtr();

    const auto lev = a_pti.GetLevel();
    const auto index = a_pti.GetPairIndex();

    m_xpold = tmp_particle_data[lev][index][TmpIdx::xold].dataPtr();
    m_ypold = tmp_particle_data[lev][index][TmpIdx::yold].dataPtr();
    m_zpold = tmp_particle_data[lev][index][TmpIdx::zold].dataPtr();
    m_uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
    m_uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
    m_uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();

    m_Phys_c = PhysConst::c;
    m_inv_c2 = 1._rt/(m_Phys_c * m_Phys_c);
}

RecordingPlaneParticleFunctor::RecordingPlaneParticleFunctor (
                        WarpXParticleContainer *pc_src, std::string species_name, int num_station_buffers)
    : m_pc_src(pc_src), m_species_name(species_name), m_num_station_buffers(num_station_buffers)
{
}

void
RecordingPlaneParticleFunctor::operator () (PinnedMemoryParticleContainer& pc_dst, int &TotalParticleCounter, int i_buffer) const
{
    amrex::Print() << " in station particle operator " << m_record_particles << "\n";
    if (!m_record_particles) return;
    auto & warpx = WarpX::GetInstance();
    auto tmp_particle_data = m_pc_src->getTmpParticleData();
    amrex::Print() << " got tmp particle ddata in operator \n";
    for (int lev = 0; lev < 1; ++lev) {
        for (WarpXParIter pti(*m_pc_src, lev); pti.isValid(); ++pti) {
            auto ptile_dst = pc_dst.DefineAndReturnParticleTile(lev, pti.index(), pti.LocalTileIndex() );
        }

        auto& particles = m_pc_src->GetParticles(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            // Temporary arrays to store copy_flag and copy_index for particles
            // that cross the recording plane
            amrex::Gpu::DeviceVector<int> FlagForPartCopy;
            amrex::Gpu::DeviceVector<int> IndexForPartCopy;

            for (WarpXParIter pti(*m_pc_src, lev); pti.isValid(); ++pti){
                auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                const auto GetParticleFilter = RecordParticles(pti, tmp_particle_data,
                                                               m_z_location);
                long const np = pti.numParticles();
                FlagForPartCopy.resize(np);
                IndexForPartCopy.resize(np);

                int* const AMREX_RESTRICT Flag = FlagForPartCopy.dataPtr();
                int* const AMREX_RESTRICT IndexLocation = IndexForPartCopy.dataPtr();

                const auto& ptile_src = particles.at(index);
                auto src_data = ptile_src.getConstParticleTileData();

                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE (int i)
                {
                    Flag[i] = GetParticleFilter(src_data, i);
                });

                const int total_partdiag_size = amrex::Scan::ExclusiveSum(np, Flag, IndexLocation);
                auto& ptile_dst = pc_dst.DefineAndReturnParticleTile(lev, pti.index(), pti.LocalTileIndex() );
                auto old_size = ptile_dst.numParticles();
                ptile_dst.resize(old_size + total_partdiag_size);
                auto dst_data = ptile_dst.getParticleTileData();
                int timecross_index = ptile_src.NumRuntimeRealComps() - 1;
                const auto GetPlaneCrossingTime = PlaneCrossingTime(pti, tmp_particle_data,
                                                  warpx.gett_new(0), m_z_location, timecross_index);
                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE (int i)
                {
                    if (Flag[i]) {
                         amrex::copyParticle(dst_data, src_data, i, old_size+IndexLocation[i]);
                         GetPlaneCrossingTime(dst_data, src_data, i, old_size+IndexLocation[i]);
                    }
                });
                amrex::Gpu::synchronize();
            }
        }
    }
    TotalParticleCounter = pc_dst.TotalNumberOfParticles();
    amrex::Print() << " end of operator : " << TotalParticleCounter << "\n";
}

void
RecordingPlaneParticleFunctor::PrepareFunctorData (const int i_buffer, bool record_particles,
                                            amrex::Real z_location, amrex::Real current_z_boost, amrex::Real tlab, int snapshot_full)
{
    amrex::ignore_unused(current_z_boost, tlab, snapshot_full);
    amrex::Print() << " record particles " << record_particles << "\n";
    m_record_particles = record_particles;
    m_z_location = z_location;
}
