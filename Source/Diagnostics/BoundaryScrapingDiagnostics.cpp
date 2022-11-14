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
#include <AMReX_ParmParse.H>

#include <set>
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
    m_varnames = {}; // No fields in boundary scraping diagnostics

    // num_buffers corresponds to the number of boundaries
    // (upper/lower domain boundary in each dimension)
    // + the EB boundary if available
    m_num_buffers = AMREX_SPACEDIM*2;
#ifdef AMREX_USE_EB
    m_num_buffers += 1;
#endif

    // Do a few checks
#ifndef WARPX_USE_OPENPMD
    amrex::Abort("You need to compile WarpX with openPMD support, in order to use BoundaryScrapingDiagnostic: -DWarpX_OPENPMD=ON");
#endif

    // Check that the output format is openPMD
    std::string error_string = std::string("You need to set `")
        .append(m_diag_name)
        .append(".format=openpmd` for the BoundaryScrapingDiagnostic.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "openpmd",
        error_string);

    // Check for the optional intervals parameter
    amrex::ParmParse pp_diag_name(m_diag_name);
    std::vector<std::string> intervals_string_vec = {"0"};
    pp_diag_name.queryarr("intervals", intervals_string_vec);
    m_intervals = utils::parser::IntervalsParser(intervals_string_vec);

}

void
BoundaryScrapingDiagnostics::InitializeFieldFunctors (int /*lev*/)
{
    // This function is usually used for field output
    // Nothing to do here for boundary scraping output,
    // since it only outputs particles
}

void
BoundaryScrapingDiagnostics::InitializeBufferData (int /*i_buffer*/, int /*lev*/)
{
    // This function is usually used for field output
    // Nothing to do here for boundary scraping output,
    // since it only outputs particles
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
            WarpXParticleContainer* pc = &mpc.GetParticleContainerFromName(species_name);
            PinnedMemoryParticleContainer* bnd_buffer = particle_buffer.getParticleBufferPointer(species_name, i_buffer);
            m_output_species[i_buffer].push_back(ParticleDiag(m_diag_name, species_name, pc, bnd_buffer));
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
BoundaryScrapingDiagnostics::DoDump (int step, int /*i_buffer*/, bool force_flush)
{
    if (force_flush) {
        return true;
    } else {
        return (m_intervals.contains(step+1));
    }
}

void
BoundaryScrapingDiagnostics::Flush (int i_buffer)
{
    auto & warpx = WarpX::GetInstance();
    ParticleBoundaryBuffer& particle_buffer = warpx.GetParticleBoundaryBuffer();

    int n_particles = 0;
    for (auto const& species_name : m_output_species_names) {
        n_particles += particle_buffer.getNumParticlesInContainer(species_name, i_buffer);
    }

    // If the saving of the particles was not set up for any of the species for this boundary
    // or if no particles have been lost, then don't write anything out.
    if (n_particles == 0) return;

    // This is not a backtransform diagnostics
    bool const isBTD = false;
    bool const isLastBTD = false;
    int const bufferID = 0;
    int const numBTDBuffers = 0;
    // The data being written out is saved in a pinned particle container
    bool const use_pinned_pc = true;
    const amrex::Geometry& geom = warpx.Geom(0); // For compatibility with `WriteToFile` ; not used

    // The data for each boundary is written out to a separate directory with the boundary name
    std::string file_prefix = m_file_prefix + "/particles_at_" + particle_buffer.boundaryName(i_buffer);

    m_flush_format->WriteToFile(
        m_varnames, m_mf_output[i_buffer], m_geom_output[i_buffer], warpx.getistep(),
        warpx.gett_new(0), m_output_species[i_buffer], nlev_output, file_prefix,
        m_file_min_digits, false, false, use_pinned_pc, isBTD,
        warpx.getistep(0), bufferID, numBTDBuffers, geom,
        isLastBTD, m_totalParticles_flushed_already[i_buffer]);

    // Now that the data has been written out, clear out the buffer
    particle_buffer.clearParticles(i_buffer);

}
