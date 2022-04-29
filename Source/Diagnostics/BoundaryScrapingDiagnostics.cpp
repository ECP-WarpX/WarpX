/* Copyright 2022 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "BoundaryScrapingDiagnostics.H"
#include "ComputeDiagFunctors/BackTransformFunctor.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "ComputeDiagFunctors/RhoFunctor.H"
#include "Diagnostics/Diagnostics.H"
#include "Diagnostics/FlushFormats/FlushFormat.H"
#include "Parallelization/WarpXCommUtil.H"
#include "ComputeDiagFunctors/BackTransformParticleFunctor.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Utils/CoarsenIO.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_BLassert.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_CoordSys.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FileSystem.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

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

    // Extract reference to the particle buffer
    ParticleBoundaryBuffer& particle_buffer = warpx.GetParticleBoundaryBuffer();

    // Initialize one ParticleDiag per species requested
    for (int i_buffer = 0; i_buffer < m_num_buffers; ++i_buffer) {
        for (auto const& species_name : m_output_species_names){
            // `particle_buffer` contains buffers for all boundaries
            // here we select the one for the EB (index: AMREX_SPACEDIM*2)
            ParticleBoundaryBuffer::BufferType& eb_buffer = particle_buffer.getParticleBuffer(species_name, AMREX_SPACEDIM*2);
            m_output_species[i_buffer].push_back(ParticleDiag(m_diag_name, species_name, nullptr, &eb_buffer));
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
    amrex::Print() << ">>> Flushing boundary scraping." << std::endl;
    auto & warpx = WarpX::GetInstance();

    m_flush_format->WriteToFile(
        m_varnames, m_mf_output[i_buffer], m_geom_output[i_buffer], warpx.getistep(),
        warpx.gett_new(0), m_output_species[i_buffer], nlev_output, m_file_prefix,
        m_file_min_digits, false, false);
}
