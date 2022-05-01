/* Copyright 2022 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "BoundaryScrapingDiagnostics.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "Diagnostics/Diagnostics.H"
#include "Diagnostics/FlushFormats/FlushFormat.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "WarpX.H"

#include <AMReX.H>

#include <string>

using namespace amrex::literals;

BoundaryScrapingDiagnostics::BoundaryScrapingDiagnostics (int i, std::string name)
    : Diagnostics(i, name)
{
    ReadParameters();
}

void
BoundaryScrapingDiagnostics::ReadParameters ()
{
    BaseReadParameters();

    // Modify some of the quantities that were initialized by default
    // in the function `BaseReadParameters`
    m_varnames_fields = {}; // No fields in boundary scraping diagnostics

    // Number of buffers = 1 for BoundaryScrapingDiagnostics.
    // (buffers are used in BTDiagnostics, and correspond to different snapshots)
    m_num_buffers = 1;
}

void
BoundaryScrapingDiagnostics::InitializeFieldFunctors (int /*lev*/)
{
}

void
BoundaryScrapingDiagnostics::InitializeBufferData (int /*i_buffer*/, int /*lev*/)
{
    // This function is usually used for field output
    // Nothing to do here for boundary scraping output, since it only outputs particles
}

void
BoundaryScrapingDiagnostics::InitializeParticleBuffer ()
{
    auto & warpx = WarpX::GetInstance();
    const MultiParticleContainer& mpc = warpx.GetPartContainer();

    // If the user does not specify any species, dump all species
    if (m_output_species_names.empty()) {
        m_output_species_names = mpc.GetSpeciesNames();
    }

    // Initialize one ParticleDiag per species requested
    ParticleBoundaryBuffer& particle_buffer = warpx.GetParticleBoundaryBuffer();
    for (int i_buffer = 0; i_buffer < m_num_buffers; ++i_buffer) {
        for (auto const& species_name : m_output_species_names){
            // `particle_buffer` contains buffers for all boundaries
            // here we select the one for the EB (index: AMREX_SPACEDIM*2)
            WarpXParticleContainer* pc = &mpc.GetParticleContainerFromName(species_name);
            PinnedMemoryParticleContainer* eb_buffer = particle_buffer.getParticleBufferPointer(species_name, AMREX_SPACEDIM*2);
            m_output_species[i_buffer].push_back(ParticleDiag(m_diag_name, species_name, pc, eb_buffer));
        }
    }
    // Initialize total number of particles flushed
    m_totalParticles_flushed_already.resize(m_num_buffers);
    for (int i_buffer = 0; i_buffer < m_num_buffers; ++i_buffer) {
        int const n_species = m_output_species_names.size();
        m_totalParticles_flushed_already[i_buffer].resize(n_species);
        for (int i_species=0; i_species<n_species; i_species++) {
            m_totalParticles_flushed_already[i_buffer][i_species] = 0;
        }
    }
}

bool
BoundaryScrapingDiagnostics::DoComputeAndPack (int /*step*/, bool /*force_flush*/)
{

    return false;
}

bool
BoundaryScrapingDiagnostics::DoDump (int /*step*/, int /*i_buffer*/, bool force_flush)
{
    if (force_flush) {
        return true;
    } else {
        return false;
    }
}

void
BoundaryScrapingDiagnostics::Flush (int i_buffer)
{
    auto & warpx = WarpX::GetInstance();

    // This is not a backtransform diagnostics, but we still set the flag `isBTD`
    // This enables:
    //   - writing the data that was accumulated in a PinnedMemoryParticleContainer
    //   - writing repeatedly to the same file
    bool const isBTD = true;
    // TODO: Change once we repeatedly call this function throughout the simulation
    bool const isLastBTD = true;
    const amrex::Geometry& geom = warpx.Geom(0); // For compatibility with `WriteToFile` ; not used

    m_flush_format->WriteToFile(
        m_varnames, m_mf_output[i_buffer], m_geom_output[i_buffer], warpx.getistep(),
        0., m_output_species[i_buffer], nlev_output, m_file_prefix,
        m_file_min_digits, false, false, isBTD, i_buffer, geom,
        isLastBTD, m_totalParticles_flushed_already[i_buffer]);
}
