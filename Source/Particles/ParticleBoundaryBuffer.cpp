/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpX.H"
#include "EmbeddedBoundary/DistanceToEB.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <ablastr/particles/NodalFieldGather.H>

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>
#include <AMReX.H>

struct IsOutsideDomainBoundary {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_phi;
    int m_idim;
    int m_iside;

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
    m_do_any_boundary.resize(numBoundaries(), 0);
    m_boundary_names.resize(numBoundaries());

    for (int i = 0; i < numBoundaries(); ++i)
    {
        m_particle_containers[i].resize(numSpecies());
        m_do_boundary_buffer[i].resize(numSpecies(), 0);
    }

#if defined(WARPX_DIM_1D_Z)
    constexpr auto idx_zlo = 0;
    constexpr auto idx_zhi = 1;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    constexpr auto idx_xlo = 0;
    constexpr auto idx_xhi = 1;
    constexpr auto idx_zlo = 2;
    constexpr auto idx_zhi = 3;
#else
    constexpr auto idx_xlo = 0;
    constexpr auto idx_xhi = 1;
    constexpr auto idx_ylo = 2;
    constexpr auto idx_yhi = 3;
    constexpr auto idx_zlo = 4;
    constexpr auto idx_zhi = 5;
#endif

    for (int ispecies = 0; ispecies < numSpecies(); ++ispecies)
    {
        amrex::ParmParse pp_species(getSpeciesNames()[ispecies]);
#if defined(WARPX_DIM_1D_Z)
        pp_species.query("save_particles_at_zlo", m_do_boundary_buffer[idx_zlo][ispecies]);
        pp_species.query("save_particles_at_zhi", m_do_boundary_buffer[idx_zhi][ispecies]);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        pp_species.query("save_particles_at_xlo", m_do_boundary_buffer[idx_xlo][ispecies]);
        pp_species.query("save_particles_at_xhi", m_do_boundary_buffer[idx_xhi][ispecies]);
        pp_species.query("save_particles_at_zlo", m_do_boundary_buffer[idx_zlo][ispecies]);
        pp_species.query("save_particles_at_zhi", m_do_boundary_buffer[idx_zhi][ispecies]);
#else
        pp_species.query("save_particles_at_xlo", m_do_boundary_buffer[idx_xlo][ispecies]);
        pp_species.query("save_particles_at_xhi", m_do_boundary_buffer[idx_xhi][ispecies]);
        pp_species.query("save_particles_at_ylo", m_do_boundary_buffer[idx_ylo][ispecies]);
        pp_species.query("save_particles_at_yhi", m_do_boundary_buffer[idx_yhi][ispecies]);
        pp_species.query("save_particles_at_zlo", m_do_boundary_buffer[idx_zlo][ispecies]);
        pp_species.query("save_particles_at_zhi", m_do_boundary_buffer[idx_zhi][ispecies]);
#endif
#ifdef AMREX_USE_EB
        pp_species.query("save_particles_at_eb", m_do_boundary_buffer[AMREX_SPACEDIM*2][ispecies]);
#endif
        // Set the flag whether the boundary is active or any species
        for (int i = 0; i < numBoundaries(); ++i) {
            if (m_do_boundary_buffer[i][ispecies]) m_do_any_boundary[i] = 1;
        }
    }

#if defined(WARPX_DIM_1D_Z)
    m_boundary_names[idx_zlo] = "zlo";
    m_boundary_names[idx_zhi] = "zhi";
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    m_boundary_names[idx_xlo] = "xlo";
    m_boundary_names[idx_xhi] = "xhi";
    m_boundary_names[idx_zlo] = "zlo";
    m_boundary_names[idx_zhi] = "zhi";
#else
    m_boundary_names[idx_xlo] = "xlo";
    m_boundary_names[idx_xhi] = "xhi";
    m_boundary_names[idx_ylo] = "ylo";
    m_boundary_names[idx_yhi] = "yhi";
    m_boundary_names[idx_zlo] = "zlo";
    m_boundary_names[idx_zhi] = "zhi";
#endif
#ifdef AMREX_USE_EB
    m_boundary_names[AMREX_SPACEDIM*2] =  "eb";
#endif

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
                amrex::Print() << Utils::TextMsg::Info(
                    "Species " + getSpeciesNames()[i] + " has "
                    + std::to_string(np) + " particles in the boundary buffer "
                    + "for side " + std::to_string(iside) + " of dim " + std::to_string(idim)
                );
            }
        }
    }
#ifdef AMREX_USE_EB
    auto& buffer = m_particle_containers[2*AMREX_SPACEDIM];
    for (int i = 0; i < numSpecies(); ++i)
    {
        int np = buffer[i].isDefined() ? buffer[i].TotalNumberOfParticles(false) : 0;
        amrex::Print() << Utils::TextMsg::Info(
            "Species " + getSpeciesNames()[i] + " has "
            + std::to_string(np) + " particles in the EB boundary buffer"
        );
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
                // do not remove particles with negative ids
                species_buffer.Redistribute(0, -1, 0, 0, false);
            }
        }
    }
}

void ParticleBoundaryBuffer::clearParticles () {
    for (int i = 0; i < numBoundaries(); ++i)
    {
        clearParticles(i);
    }
}

void ParticleBoundaryBuffer::clearParticles (int const i) {
    auto& buffer = m_particle_containers[i];
    for (int ispecies = 0; ispecies < numSpecies(); ++ispecies)
    {
        auto& species_buffer = buffer[ispecies];
        if (species_buffer.isDefined()) species_buffer.clearParticles();
    }
}

void ParticleBoundaryBuffer::gatherParticles (MultiParticleContainer& mypc,
                                              const amrex::Vector<const amrex::MultiFab*>& distance_to_eb)
{
    WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles");

    using PIter = amrex::ParConstIter<0,0,PIdx::nattribs>;
    const auto& warpx_instance = WarpX::GetInstance();
    const amrex::Geometry& geom = warpx_instance.Geom(0);
    auto plo = geom.ProbLoArray();
    auto phi = geom.ProbHiArray();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (geom.isPeriodic(idim)) continue;
        for (int iside = 0; iside < 2; ++iside)
        {
            auto& buffer = m_particle_containers[2*idim+iside];
            for (int i = 0; i < numSpecies(); ++i)
            {
                if (!m_do_boundary_buffer[2*idim+iside][i]) continue;
                const WarpXParticleContainer& pc = mypc.GetParticleContainer(i);
                if (!buffer[i].isDefined())
                {
                    buffer[i] = pc.make_alike<amrex::PinnedArenaAllocator>();
                    buffer[i].AddIntComp("timestamp", false);
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

                        auto predicate = IsOutsideDomainBoundary{plo, phi, idim, iside};

                        const auto ptile_data = ptile.getConstParticleTileData();

                        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
                        amrex::ReduceData<int> reduce_data(reduce_op);
                        {
                          WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles::count_out_of_bounds");
                          amrex::RandomEngine rng{};
                          reduce_op.eval(np, reduce_data, [=] AMREX_GPU_HOST_DEVICE (int ip)
                                         { return predicate(ptile_data, ip, rng) ? 1 : 0; });
                        }

                        auto dst_index = ptile_buffer.numParticles();
                        {
                          WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles::resize");
                          ptile_buffer.resize(dst_index + amrex::get<0>(reduce_data.value()));
                        }
                        {
                          WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles::filterAndTransform");
                          int timestamp_index = ptile_buffer.NumRuntimeIntComps()-1;
                          int timestep = warpx_instance.getistep(0);

                          amrex::filterAndTransformParticles(ptile_buffer, ptile,
                                                             predicate,
                                                             CopyAndTimestamp{timestamp_index, timestep},
                                                             0, dst_index);
                        }
                    }
                }
            }
        }
    }

#ifdef AMREX_USE_EB
    WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles::EB");

    auto& buffer = m_particle_containers[m_particle_containers.size()-1];
    for (int i = 0; i < numSpecies(); ++i)
    {
        if (!m_do_boundary_buffer[AMREX_SPACEDIM*2][i]) continue;
        const auto& pc = mypc.GetParticleContainer(i);
        if (!buffer[i].isDefined())
        {
            buffer[i] = pc.make_alike<amrex::PinnedArenaAllocator>();
            buffer[i].AddIntComp("timestamp", false);
        }
        auto& species_buffer = buffer[i];
        for (int lev = 0; lev < pc.numLevels(); ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            auto dxi = warpx_instance.Geom(lev).InvCellSizeArray();
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

                using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
                auto predicate = [=] AMREX_GPU_HOST_DEVICE (const SrcData& /*src*/, const int ip)
                /* NVCC 11.3.109 chokes in C++17 on this: noexcept */
                  {
                    amrex::ParticleReal xp, yp, zp;
                    getPosition(ip, xp, yp, zp);

                    amrex::Real phi_value  = ablastr::particles::doGatherScalarFieldNodal(
                                                                      xp, yp, zp, phiarr, dxi, plo
                                                                      );
                    return phi_value < 0.0 ? 1 : 0;
                  };

                const auto ptile_data = ptile.getConstParticleTileData();

                amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
                amrex::ReduceData<int> reduce_data(reduce_op);
                {
                  WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles::count_out_of_boundsEB");
                  reduce_op.eval(np, reduce_data, [=] AMREX_GPU_HOST_DEVICE (int ip)
                                 { return predicate(ptile_data, ip) ? 1 : 0; });
                }

                auto dst_index = ptile_buffer.numParticles();
                {
                  WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles::resize_eb");
                  ptile_buffer.resize(dst_index + amrex::get<0>(reduce_data.value()));
                }

                int timestamp_index = ptile_buffer.NumRuntimeIntComps()-1;
                int timestep = warpx_instance.getistep(0);
                {
                  WARPX_PROFILE("ParticleBoundaryBuffer::gatherParticles::filterTransformEB");
                  amrex::filterAndTransformParticles(ptile_buffer, ptile, predicate,
                                                     CopyAndTimestamp{timestamp_index, timestep}, 0, dst_index);
                }
            }
        }
    }
#else
    amrex::ignore_unused(distance_to_eb);
#endif
}

int ParticleBoundaryBuffer::getNumParticlesInContainer(
        const std::string species_name, int boundary, bool local) {

    auto& buffer = m_particle_containers[boundary];
    auto index = WarpX::GetInstance().GetPartContainer().getSpeciesID(species_name);

    if (buffer[index].isDefined()) return buffer[index].TotalNumberOfParticles(false, local);
    else return 0;
}

PinnedMemoryParticleContainer &
ParticleBoundaryBuffer::getParticleBuffer(const std::string species_name, int boundary) {

    auto& buffer = m_particle_containers[boundary];
    auto index = WarpX::GetInstance().GetPartContainer().getSpeciesID(species_name);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_do_boundary_buffer[boundary][index],
                                     "Attempted to get particle buffer for boundary "
                                     + std::to_string(boundary) + ", which is not used!");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(buffer[index].isDefined(),
                                     "Tried to get a buffer that is not defined!");

    return buffer[index];
}

PinnedMemoryParticleContainer *
ParticleBoundaryBuffer::getParticleBufferPointer(const std::string species_name, int boundary) {

    auto& buffer = m_particle_containers[boundary];
    auto index = WarpX::GetInstance().GetPartContainer().getSpeciesID(species_name);

    return &buffer[index];
}
