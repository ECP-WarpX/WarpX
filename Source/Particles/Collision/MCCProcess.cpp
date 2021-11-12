/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MCCProcess.H"
#include "WarpX.H"

MCCProcess::MCCProcess (
                        const std::string& scattering_process,
                        const std::string& cross_section_file,
                        const amrex::Real energy )
{
    amrex::Print() << "Reading file " << cross_section_file << " for "
                   << scattering_process << " scattering cross-sections.\n";

    // read the cross-section data file into memory
    readCrossSectionFile(cross_section_file, m_energies, m_sigmas_h);

    init(scattering_process, energy);
}

template <typename InputVector>
MCCProcess::MCCProcess (
                        const std::string& scattering_process,
                        const InputVector&& energies,
                        const InputVector&& sigmas,
                        const amrex::Real energy )
{
    m_energies.insert(m_energies.begin(), std::begin(energies), std::end(energies));
    m_sigmas_h.insert(m_sigmas_h.begin(), std::begin(sigmas),   std::end(sigmas));

    init(scattering_process, energy);
}

void
MCCProcess::init (const std::string& scattering_process, const amrex::Real energy)
{
    m_exe_h.m_sigmas_data = m_sigmas_h.data();

    // save energy grid parameters for easy use
    m_grid_size = m_energies.size();
    m_exe_h.m_energy_lo = m_energies[0];
    m_exe_h.m_energy_hi = m_energies[m_grid_size-1];
    m_exe_h.m_sigma_lo = m_sigmas_h[0];
    m_exe_h.m_sigma_hi = m_sigmas_h[m_grid_size-1];
    m_exe_h.m_dE = (m_exe_h.m_energy_hi - m_exe_h.m_energy_lo)/(m_grid_size - 1.);
    m_exe_h.m_energy_penalty = energy;
    m_exe_h.m_type = parseProcessType(scattering_process);

    // sanity check cross-section energy grid
    sanityCheckEnergyGrid(m_energies, m_exe_h.m_dE);

    // check that the cross-section is 0 at the energy cost if the energy
    // cost is > 0 - this is to prevent the possibility of negative left
    // over energy after a collision event
    if (m_exe_h.m_energy_penalty > 0) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            (getCrossSection(m_exe_h.m_energy_penalty) == 0),
            "Cross-section > 0 at energy cost for collision."
        );
    }

#ifdef AMREX_USE_GPU
    m_exe_d = m_exe_h;
    m_sigmas_d.resize(m_sigmas_h.size());
    m_exe_d.m_sigmas_data = m_sigmas_d.data();
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_sigmas_h.begin(), m_sigmas_h.end(),
                          m_sigmas_d.begin());
    amrex::Gpu::streamSynchronize();
#endif
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
    } else if (scattering_process == "ionization") {
        return MCCProcessType::IONIZATION;
    } else if (scattering_process.find("excitation") != std::string::npos) {
        return MCCProcessType::EXCITATION;
    } else {
        return MCCProcessType::INVALID;
    }
}

void
MCCProcess::readCrossSectionFile (
                                  const std::string cross_section_file,
                                  amrex::Vector<amrex::Real>& energies,
                                  amrex::Gpu::HostVector<amrex::Real>& sigmas )
{
    std::ifstream infile(cross_section_file);
    if(!infile.is_open()) amrex::Abort("Failed to open cross-section data file");

    double energy, sigma;
    while (infile >> energy >> sigma) {
        energies.push_back(energy);
        sigmas.push_back(sigma);
    }
    if (infile.bad()) amrex::Abort("Failed to read cross-section data from file.");
    infile.close();
}

void
MCCProcess::sanityCheckEnergyGrid (
                                   const amrex::Vector<amrex::Real>& energies,
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
