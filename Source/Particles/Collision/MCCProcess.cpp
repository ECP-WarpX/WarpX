/* Copyright 2021 Roelof Groenewald
 *
 * This file is part of WarpX.
 *
 * License: ????
 */
#include "MCCProcess.H"
#include "WarpX.H"

MCCProcess::MCCProcess (
                        const std::string& scattering_process,
                        const std::string& cross_section_file,
                        const amrex::Real energy )
    : type(parseProcessType(scattering_process))
{
    amrex::Print() << "Reading file " << cross_section_file << " for "
                   << scattering_process << " scattering cross-sections.\n";

    // read the cross-section data file into memory
    readCrossSectionFile(cross_section_file, m_energies, m_sigmas);
    m_energies_data = m_energies.data();
    m_sigmas_data = m_sigmas.data();

    // save energy grid parameters for easy use
    m_grid_size = m_energies.size();
    m_energy_lo = m_energies[0];
    m_energy_hi = m_energies[m_grid_size-1];
    m_sigma_lo = m_sigmas[0];
    m_sigma_hi = m_sigmas[m_grid_size-1];
    m_dE = m_energies[1] - m_energies[0];

    energy_penalty = energy;

    // sanity check cross-section energy grid
    sanityCheckEnergyGrid(m_energies, m_dE);
}

MCCProcessType
MCCProcess::parseProcessType(const std::string& scattering_process)
{
    if (scattering_process == "elastic") {
        return MCCProcessType::ELASTIC;
    } else if (scattering_process == "back") {
        return MCCProcessType::BACK;
    } else if (scattering_process == "charge_exchange") {
        return MCCProcessType::CHARGE_EXCHANGE;
    } else if (scattering_process.find("excitation") != std::string::npos) {
        return MCCProcessType::EXCITATION;
    } else {
        return MCCProcessType::INVALID;
    }
}


void
MCCProcess::readCrossSectionFile (
                                  const std::string cross_section_file,
                                  MCCProcess::VectorType<amrex::Real>& energies,
                                  MCCProcess::VectorType<amrex::Real>& sigmas )
{
    std::ifstream infile(cross_section_file);
    double energy, sigma;
    while (infile >> energy >> sigma) {
        energies.push_back(energy);
        sigmas.push_back(sigma);
    }
}

void
MCCProcess::sanityCheckEnergyGrid (
                                   const VectorType<amrex::Real>& energies,
                                   amrex::Real dE
                                   )
{
    // confirm that the input data for the cross-section was provided with
    // equal energy steps, otherwise the linear interpolation will fail
    for (unsigned i = 1; i < energies.size(); i++) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                                         (std::abs(energies[i] - energies[i-1] - dE) < dE / 100.0),
                                         "Energy grid not evenly spaced."
                                         );
    }
}
