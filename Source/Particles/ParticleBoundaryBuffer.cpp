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
    int num_boundaries = AMREX_SPACEDIM*2;  // two for each dim, hi and lo
#ifdef AMREX_USE_EB
    num_boundaries += 1;  // only one of these for now
#endif
    m_particle_containers.resize(num_boundaries);

    amrex::ParmParse pp_particles("particles");
    pp_particles.queryarr("species_names", m_species_names);
    const int nspecies = m_species_names.size();

    for (int i = 0; i < num_boundaries; ++i)
    {
        m_particle_containers[i].resize(nspecies);
    }

    m_do_boundary_buffer.resize(num_boundaries);
    for (int i = 0; i < num_boundaries; ++i)
    {
        m_do_boundary_buffer[i].resize(nspecies, 0);
    }

    for (int i = 0; i < num_boundaries; ++i)
    {
        for (int ispecies = 0; ispecies < nspecies; ++ispecies)
        {
            amrex::ParmParse pp_species(m_species_names[ispecies]);
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
    const int nspecies = nSpecies();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        for (int iside = 0; iside < 2; ++iside)
        {
            auto& buffer = m_particle_containers[2*idim+iside];
            for (int i = 0; i < nspecies; ++i)
            {
                int np = buffer[i].isDefined() ? buffer[i].TotalNumberOfParticles(false) : 0;
                amrex::Print() << "Species " << m_species_names[i] << " has "
                               << np << " particles in the boundary buffer "
                               << " for side " << iside << " of dim " << idim << "\n";
            }
        }
    }
#ifdef AMREX_USE_EB
    auto& buffer = m_particle_containers[2*AMREX_SPACEDIM];
    for (int i = 0; i < nspecies; ++i)
    {
        int np = buffer[i].isDefined() ? buffer[i].TotalNumberOfParticles(false) : 0;
        amrex::Print() << "Species " << m_species_names[i] << " has "
                       << np << " particles in the EB boundary buffer \n";
    }
#endif
}

void ParticleBoundaryBuffer::gatherParticles (MultiParticleContainer& mypc,
                                              const amrex::Vector<const amrex::MultiFab*>& distance_to_eb) {
    const auto& warpx_instance = WarpX::GetInstance();
    const amrex::Geometry& geom = warpx_instance.Geom(0);
    const int nspecies = mypc.nSpecies();
    auto plo = geom.ProbLoArray();
    auto phi = geom.ProbHiArray();
    auto dxi = geom.InvCellSizeArray();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (geom.isPeriodic(idim)) continue;
        for (int iside = 0; iside < 2; ++iside)
        {
            auto& buffer = m_particle_containers[2*idim+iside];
            for (int i = 0; i < nspecies; ++i)
            {
                if (!m_do_boundary_buffer[2*idim+iside][i]) continue;
                const auto& pc = mypc.GetParticleContainer(i);
                if (!buffer[i].isDefined())
                {
                    buffer[i] = ParticleBuffer::getTmpPC<amrex::PinnedArenaAllocator>(&pc);
                }
                auto& species_buffer = buffer[i];
                species_buffer.addParticles(pc, IsOutsideDomainBoundary{plo, phi, idim, iside}, true);
            }
        }
    }

#ifdef AMREX_USE_EB
    auto& buffer = m_particle_containers[m_particle_containers.size()-1];
    for (int i = 0; i < nspecies; ++i)
    {
        const auto& pc = mypc.GetParticleContainer(i);
        if (!buffer[i].isDefined())
        {
            buffer[i] = ParticleBuffer::getTmpPC<amrex::PinnedArenaAllocator>(&pc);
        }
        auto& species_buffer = buffer[i];
        for (int lev = 0; lev < pc.numLevels(); ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(amrex::ParConstIter<0,0,PIdx::nattribs> pti(pc, lev); pti.isValid(); ++pti)
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

                using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
                auto count = amrex::filterParticles(ptile_buffer, ptile,
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
                                                    0, dst_index, np);
                ptile_buffer.resize(dst_index + count);
            }
        }
    }
#else
    amrex::ignore_unused(distance_to_eb, dxi);
#endif
}
