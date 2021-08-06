/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */


#include "EmbeddedBoundary/DistanceToEB.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

ParticleBoundaryBuffer::ParticleBoundaryBuffer ()
{
    m_particle_containers.resize(numBoundaries());

    for (int i = 0; i < numBoundaries(); ++i)
    {
        m_particle_containers[i].resize(numSpecies());
    }

    m_do_boundary_buffer.resize(numBoundaries());
    for (int i = 0; i < numBoundaries(); ++i)
    {
        m_do_boundary_buffer[i].resize(numSpecies(), 0);
    }

    for (int i = 0; i < numBoundaries(); ++i)
    {
        for (int ispecies = 0; ispecies < numSpecies(); ++ispecies)
        {
            amrex::ParmParse pp_species(getSpeciesNames()[ispecies]);
            pp_species.query("scrape_xlo", m_do_boundary_buffer[0][ispecies]);
            pp_species.query("scrape_xhi", m_do_boundary_buffer[1][ispecies]);
#if AMREX_SPACEDIM == 2
            pp_species.query("scrape_zlo", m_do_boundary_buffer[2][ispecies]);
            pp_species.query("scrape_zhi", m_do_boundary_buffer[3][ispecies]);
#else
            pp_species.query("scrape_ylo", m_do_boundary_buffer[2][ispecies]);
            pp_species.query("scrape_yhi", m_do_boundary_buffer[3][ispecies]);
            pp_species.query("scrape_zlo", m_do_boundary_buffer[4][ispecies]);
            pp_species.query("scrape_zhi", m_do_boundary_buffer[5][ispecies]);
#endif
#ifdef AMREX_USE_EB
            pp_species.query("scrape_eb", m_do_boundary_buffer[AMREX_SPACEDIM*2][ispecies]);
#endif
        }
    }
}

void ParticleBoundaryBuffer::printNumParticles () const {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        for (int iside = 0; iside < 2; ++iside)
        {
            auto& buffer = m_particle_containers[2*idim+iside];
            for (int i = 0; i < numSpecies(); ++i)
            {
                int np = buffer[i].isDefined() ? buffer[i].TotalNumberOfParticles(false) : 0;
                amrex::Print() << "Species " << getSpeciesNames()[i] << " has "
                               << np << " particles in the boundary buffer "
                               << " for side " << iside << " of dim " << idim << "\n";
            }
        }
    }
#ifdef AMREX_USE_EB
    auto& buffer = m_particle_containers[2*AMREX_SPACEDIM];
    for (int i = 0; i < numSpecies(); ++i)
    {
        int np = buffer[i].isDefined() ? buffer[i].TotalNumberOfParticles(false) : 0;
        amrex::Print() << "Species " << getSpeciesNames()[i] << " has "
                       << np << " particles in the EB boundary buffer \n";
    }
#endif
}

void ParticleBoundaryBuffer::clearParticles () {
    for (int i = 0; i < numBoundaries(); ++i)
    {
        auto& buffer = m_particle_containers[i];
        if (buffer[i].isDefined())
        {
            buffer[i].clearParticles();
        }
    }
}

void ParticleBoundaryBuffer::gatherParticles (MultiParticleContainer& mypc,
                                              const amrex::Vector<const amrex::MultiFab*>& distance_to_eb) {
    using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
    using DstData = ParticleBuffer::BufferType<amrex::PinnedArenaAllocator>::ParticleTileType::ParticleTileDataType;
    using PIter = amrex::ParConstIter<0,0,PIdx::nattribs>;
    const auto& warpx_instance = WarpX::GetInstance();
    const amrex::Geometry& geom = warpx_instance.Geom(0);
    auto plo = geom.ProbLoArray();
    auto phi = geom.ProbHiArray();
    auto dxi = geom.InvCellSizeArray();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (geom.isPeriodic(idim)) continue;
        for (int iside = 0; iside < 2; ++iside)
        {
            auto& buffer = m_particle_containers[2*idim+iside];
            for (int i = 0; i < numSpecies(); ++i)
            {
                if (!m_do_boundary_buffer[2*idim+iside][i]) continue;
                const auto& pc = mypc.GetParticleContainer(i);
                if (!buffer[i].isDefined())
                {
                    buffer[i] = ParticleBuffer::getTmpPC<amrex::PinnedArenaAllocator>(&pc);
                    buffer[i].AddRealComp(false);  // for timestamp
                }
                auto& species_buffer = buffer[i];
                for (int lev = 0; lev < pc.numLevels(); ++lev)
                {
                    const auto& plevel = pc.GetParticles(lev);
                    for(PIter pti(pc, lev); pti.isValid(); ++pti)
                    {
                        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                        if(plevel.find(index) == plevel.end()) continue;

                        auto& ptile_buffer = species_buffer.DefineAndReturnParticleTile(
                                                        lev, pti.index(), pti.LocalTileIndex());
                        const auto& ptile = plevel.at(index);
                        auto np = ptile.numParticles();
                        if (np == 0) continue;

                        auto dst_index = ptile_buffer.numParticles();
                        ptile_buffer.resize(dst_index + np);

                        int timestamp_index = ptile_buffer.NumRuntimeRealComps()-1;
                        amrex::Real time = warpx_instance.gett_new(0);
                        auto count = amrex::filterAndTransformParticles(ptile_buffer, ptile,
                            IsOutsideDomainBoundary{plo, phi, idim, iside},
                            [=] AMREX_GPU_HOST_DEVICE (const DstData& dst, const SrcData& src,
                                                       int src_i, int dst_i) noexcept
                        {
                            dst.m_aos[dst_i] = src.m_aos[src_i];
                            for (int j = 0; j < SrcData::NAR; ++j)
                                dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
                            for (int j = 0; j < src.m_num_runtime_real; ++j)
                                dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
                            for (int j = 0; j < src.m_num_runtime_int; ++j)
                                dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
                            dst.m_runtime_rdata[timestamp_index][dst_i] = time;
                        },
                                                                        0, dst_index);
                        ptile_buffer.resize(dst_index + count);
                    }
                }
            }
        }
    }

#ifdef AMREX_USE_EB
    auto& buffer = m_particle_containers[m_particle_containers.size()-1];
    for (int i = 0; i < numSpecies(); ++i)
    {
        const auto& pc = mypc.GetParticleContainer(i);
        if (!buffer[i].isDefined())
        {
            buffer[i] = ParticleBuffer::getTmpPC<amrex::PinnedArenaAllocator>(&pc);
            buffer[i].AddRealComp(false);  // for timestamp
        }
        auto& species_buffer = buffer[i];
        for (int lev = 0; lev < pc.numLevels(); ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(PIter pti(pc, lev); pti.isValid(); ++pti)
            {
                auto phiarr = (*distance_to_eb[lev])[pti].array();  // signed distance function
                auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                if(plevel.find(index) == plevel.end()) continue;

                const auto getPosition = GetParticlePosition(pti);
                auto& ptile_buffer = species_buffer.DefineAndReturnParticleTile(lev, pti.index(),
                                                                                pti.LocalTileIndex());
                const auto& ptile = plevel.at(index);
                auto np = ptile.numParticles();
                if (np == 0) continue;

                auto dst_index = ptile_buffer.numParticles();
                ptile_buffer.resize(dst_index + np);

                int timestamp_index = ptile_buffer.NumRuntimeRealComps()-1;
                amrex::Real time = warpx_instance.gett_new(0);
                auto count = amrex::filterAndTransformParticles(ptile_buffer, ptile,
                    [=] AMREX_GPU_HOST_DEVICE (const SrcData& /*src*/, const int ip) noexcept
                    {
                        amrex::ParticleReal xp, yp, zp;
                        getPosition(ip, xp, yp, zp);

                        int ii, jj, kk;
                        amrex::Real W[AMREX_SPACEDIM][2];
                        DistanceToEB::compute_weights(xp, yp, zp, plo, dxi, ii, jj, kk, W);

                        amrex::Real phi_value  = DistanceToEB::interp_distance(ii, jj, kk, W, phiarr);
                        return phi_value < 0.0 ? 1 : 0;
                    },
                            [=] AMREX_GPU_HOST_DEVICE (const DstData& dst, const SrcData& src,
                                                       int src_i, int dst_i) noexcept
                        {
                            dst.m_aos[dst_i] = src.m_aos[src_i];
                            for (int j = 0; j < SrcData::NAR; ++j)
                                dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
                            for (int j = 0; j < src.m_num_runtime_real; ++j)
                                dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
                            for (int j = 0; j < src.m_num_runtime_int; ++j)
                                dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
                            dst.m_runtime_rdata[timestamp_index][dst_i] = time;
                        },
                                                    0, dst_index);
                ptile_buffer.resize(dst_index + count);
            }
        }
    }
#else
    amrex::ignore_unused(distance_to_eb, dxi);
#endif
}
