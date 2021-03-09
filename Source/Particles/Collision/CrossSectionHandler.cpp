/* Copyright 2021 Roelof Groenewald
 *
 * This file is part of WarpX.
 *
 * License: ????
 */
#include "CrossSectionHandler.H"
#include "WarpX.H"

// using namespace amrex::literals;

CrossSectionHandler::CrossSectionHandler (
        const std::string scattering_process,
        const std::string cross_section_file )
{
    amrex::Print() << "Reading file " << cross_section_file << " for the "
        << scattering_process << " scattering cross-sections.\n";

    name = scattering_process;

    // read the cross-section data file into memory
    readCrossSectionFile(cross_section_file, m_energies, m_sigmas);

    // sanity check cross-section energy grid
    // sanityCheckEnergyGrid(m_energies);

    // save energy grid parameters for easy use
    m_grid_size = m_energies.size();
    m_energy_lo = m_energies[0];
    m_energy_hi = m_energies[m_grid_size-1];
    m_sigma_lo = m_sigmas[0];
    m_sigma_hi = m_sigmas[m_grid_size-1];
    m_dE = m_energies[1] - m_energies[0];

    // #TODO should fix this
    energy_penalty = 0.0;
}

amrex::Real
CrossSectionHandler::getCrossSection ( amrex::Real E_coll ) const
{
    if (E_coll < m_energy_lo)
    {
        return m_sigma_lo;
    }
    else if (E_coll > m_energy_hi)
    {
        return m_sigma_hi;
    }
    else
    {
        // calculate index of bounding energy pairs
        amrex::Real temp = (E_coll - m_energy_lo) / m_dE;
        int idx_1 = std::floor(temp);
        int idx_2 = std::ceil(temp);

        // linearly interpolate to the given energy value
        temp -= idx_1;
        return (
            m_sigmas[idx_1] + (m_sigmas[idx_2] - m_sigmas[idx_1]) * temp
        );
    }
}

void
CrossSectionHandler::readCrossSectionFile (
    const std::string cross_section_file,
    amrex::Vector<amrex::Real>& energies,
    amrex::Vector<amrex::Real>& sigmas )
{
    // as a test we'll just use fixed values
    energies.resize(3);
    sigmas.resize(3);
    energies[0] = 0.0;
    energies[1] = 0.5;
    energies[2] = 1.0;
    sigmas[0] = 5e-18;
    sigmas[1] = 2e-18;
    sigmas[2] = 1e-18;
}

void
CrossSectionHandler::sanityCheckEnergyGrid (
    const amrex::Vector<amrex::Real> energies )
{
    // TODO insert sanity checking logic specifically checking
    // for constant dE between points
    amrex::Print() << energies[0] << " sanity checking data\n";
}
