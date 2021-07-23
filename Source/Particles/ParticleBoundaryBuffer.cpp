/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

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
		for (int j = 0; j < nspecies; ++j)
		{
			m_particle_containers[i][j].isDefined();
		}
	}

    m_do_boundary_buffer.resize(AMREX_SPACEDIM*2);
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
        }
    }
}

void ParticleBoundaryBuffer::printNumParticles () const {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        for (int iside = 0; iside < 2; ++iside)
        {
            auto& buffer = m_particle_containers[2*idim+iside];
            const int nspecies = nSpecies();
            for (int i = 0; i < nspecies; ++i)
            {
                int np;
                if (buffer[i].isDefined()) { np = buffer[i].TotalNumberOfParticles(false); }
                else { np = 0; }
                amrex::Print() << "Species " << m_species_names[i] << " has "
                               << np << " particles in the boundary buffer "
                               << " for side " << iside << " of dim " << idim << "\n";
            }
        }
    }
}

void ParticleBoundaryBuffer::gatherParticles (MultiParticleContainer& mypc) {
    const auto& warpx_instance = WarpX::GetInstance();
    const amrex::Geometry& geom = warpx_instance.Geom(0);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (geom.isPeriodic(idim)) continue;
        for (int iside = 0; iside < 2; ++iside)
        {
            auto plo = geom.ProbLoArray();
            auto phi = geom.ProbHiArray();
            auto& buffer = m_particle_containers[2*idim+iside];
            const int nspecies = mypc.nSpecies();
            for (int i = 0; i < nspecies; ++i)
            {
                if (!m_do_boundary_buffer[2*idim+iside][i]) continue;
                const auto& pc = mypc.GetParticleContainer(i);
                if (!buffer[i].isDefined())
                {
                    buffer[i] = ParticleBuffer::getTmpPC<amrex::PinnedArenaAllocator>(&pc);
                }
                auto& species_buffer = buffer[i];
                using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
                species_buffer.addParticles(pc,
                                            [=] AMREX_GPU_DEVICE (const SrcData& src, int ip, const amrex::RandomEngine& /*engine*/)
                {
                    const auto& p = src.getSuperParticle(ip);
                    if (iside == 0) {
                        if (p.pos(idim) < plo[idim]) { return 1; }
                    } else {
                        if (p.pos(idim) >= phi[idim]) { return 1; }
                    }
                    return 0;
                }, true);
            }
        }
    }
}
