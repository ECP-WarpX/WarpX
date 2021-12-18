/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */


#include "EmbeddedBoundary/DistanceToEB.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/Gather/ScalarFieldGather.H"

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

struct IsOutsideDomainBoundary {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_phi;
    int m_idim;
    int m_iside;

    template <typename SrcData>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int operator() (const SrcData& src,
                    int ip, const amrex::RandomEngine& /*engine*/) const noexcept
    {
        const auto& p = src.getSuperParticle(ip);
        if (m_iside == 0) {
            if (p.pos(m_idim) < m_plo[m_idim]) { return 1; }
        } else {
            if (p.pos(m_idim) >= m_phi[m_idim]) { return 1; }
        }
        return 0;
    }
};

struct CopyAndTimestamp {
    int m_index;
    int m_step;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator() (const DstData& dst, const SrcData& src,
                     int src_i, int dst_i) const noexcept
    {
        dst.m_aos[dst_i] = src.m_aos[src_i];
        for (int j = 0; j < SrcData::NAR; ++j)
            dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
        for (int j = 0; j < src.m_num_runtime_real; ++j)
            dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
        for (int j = 0; j < src.m_num_runtime_int; ++j)
            dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
        dst.m_runtime_idata[m_index][dst_i] = m_step;
    }
};

ParticleBoundaryBuffer::ParticleBoundaryBuffer ()
{
    m_particle_containers.resize(numBoundaries());
    m_do_boundary_buffer.resize(numBoundaries());

    for (int i = 0; i < numBoundaries(); ++i)
    {
        m_particle_containers[i].resize(numSpecies());
        m_do_boundary_buffer[i].resize(numSpecies(), 0);
    }

    for (int ispecies = 0; ispecies < numSpecies(); ++ispecies)
    {
        amrex::ParmParse pp_species(getSpeciesNames()[ispecies]);
#if defined(WARPX_DIM_1D_Z)
        pp_species.query("save_particles_at_zlo", m_do_boundary_buffer[0][ispecies]);
        pp_species.query("save_particles_at_zhi", m_do_boundary_buffer[1][ispecies]);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        pp_species.query("save_particles_at_xlo", m_do_boundary_buffer[0][ispecies]);
        pp_species.query("save_particles_at_xhi", m_do_boundary_buffer[1][ispecies]);
        pp_species.query("save_particles_at_zlo", m_do_boundary_buffer[2][ispecies]);
        pp_species.query("save_particles_at_zhi", m_do_boundary_buffer[3][ispecies]);
#else
        pp_species.query("save_particles_at_xlo", m_do_boundary_buffer[0][ispecies]);
        pp_species.query("save_particles_at_xhi", m_do_boundary_buffer[1][ispecies]);
        pp_species.query("save_particles_at_ylo", m_do_boundary_buffer[2][ispecies]);
        pp_species.query("save_particles_at_yhi", m_do_boundary_buffer[3][ispecies]);
        pp_species.query("save_particles_at_zlo", m_do_boundary_buffer[4][ispecies]);
        pp_species.query("save_particles_at_zhi", m_do_boundary_buffer[5][ispecies]);
#endif
#ifdef AMREX_USE_EB
        pp_species.query("save_particles_at_eb", m_do_boundary_buffer[AMREX_SPACEDIM*2][ispecies]);
#endif
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
                               << "for side " << iside << " of dim " << idim << "\n";
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

void ParticleBoundaryBuffer::redistribute () {
    for (int i = 0; i < numBoundaries(); ++i)
    {
        auto& buffer = m_particle_containers[i];
        for (int ispecies = 0; ispecies < numSpecies(); ++ispecies)
        {
            auto& species_buffer = buffer[ispecies];
            if (species_buffer.isDefined()) {
                species_buffer.Redistribute();
            }
        }
    }
}

void ParticleBoundaryBuffer::clearParticles () {
    for (int i = 0; i < numBoundaries(); ++i)
    {
        auto& buffer = m_particle_containers[i];
        for (int ispecies = 0; ispecies < numSpecies(); ++ispecies)
        {
            auto& species_buffer = buffer[ispecies];
            if (species_buffer.isDefined()) species_buffer.clearParticles();
        }
    }
}

void ParticleBoundaryBuffer::gatherParticles (MultiParticleContainer& mypc,
                                              const amrex::Vector<const amrex::MultiFab*>& distance_to_eb)
{
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
                    buffer[i].AddIntComp(false);  // for timestamp
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

                        int timestamp_index = ptile_buffer.NumRuntimeIntComps()-1;
                        int timestep = warpx_instance.getistep(0);
                        auto count = amrex::filterAndTransformParticles(ptile_buffer, ptile,
                                             IsOutsideDomainBoundary{plo, phi, idim, iside},
                                                CopyAndTimestamp{timestamp_index, timestep},
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
            buffer[i].AddIntComp(false);  // for timestamp
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

                int timestamp_index = ptile_buffer.NumRuntimeIntComps()-1;
                int timestep = warpx_instance.getistep(0);
                using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
                auto count = amrex::filterAndTransformParticles(ptile_buffer, ptile,
                    [=] AMREX_GPU_HOST_DEVICE (const SrcData& /*src*/, const int ip)
                    /* NVCC 11.3.109 chokes in C++17 on this: noexcept */
                    {
                        amrex::ParticleReal xp, yp, zp;
                        getPosition(ip, xp, yp, zp);

                        amrex::Real phi_value  = doGatherScalarFieldNodal(
                            xp, yp, zp, phiarr, dxi, plo
                        );
                        return phi_value < 0.0 ? 1 : 0;
                    },
                    CopyAndTimestamp{timestamp_index, timestep}, 0, dst_index);
                ptile_buffer.resize(dst_index + count);
            }
        }
    }
#else
    amrex::ignore_unused(distance_to_eb, dxi);
#endif
}

int ParticleBoundaryBuffer::getNumParticlesInContainer(
        const std::string species_name, int boundary) {

    auto& buffer = m_particle_containers[boundary];
    auto index = WarpX::GetInstance().GetPartContainer().getSpeciesID(species_name);

    if (buffer[index].isDefined()) return buffer[index].TotalNumberOfParticles(false);
    else return 0;
}

ParticleBuffer::BufferType<amrex::PinnedArenaAllocator>&
ParticleBoundaryBuffer::getParticleBuffer(const std::string species_name, int boundary) {

    auto& buffer = m_particle_containers[boundary];
    auto index = WarpX::GetInstance().GetPartContainer().getSpeciesID(species_name);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_do_boundary_buffer[boundary][index],
                                     "Attempted to get particle buffer for boundary "
                                     + boundary + ", which is not used!");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(buffer[index].isDefined(),
                                     "Tried to get a buffer that is not defined!");

    return buffer[index];
}
