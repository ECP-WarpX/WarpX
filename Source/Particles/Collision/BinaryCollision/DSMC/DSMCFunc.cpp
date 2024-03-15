/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#include "DSMCFunc.H"

/**
 * \brief Constructor of the DSMCFunc class
 *
 * @param[in] collision_name the name of the collision
 * @param[in] mypc pointer to the MultiParticleContainer
 * @param[in] isSameSpecies whether the two colliding species are the same
 */
DSMCFunc::DSMCFunc (
    const std::string collision_name,
    [[maybe_unused]] MultiParticleContainer const * const mypc,
    [[maybe_unused]] const bool isSameSpecies )
{
    using namespace amrex::literals;

    const amrex::ParmParse pp_collision_name(collision_name);

    // query for a list of collision processes
    // these could be elastic, excitation, charge_exchange, back, etc.
    amrex::Vector<std::string> scattering_process_names;
    pp_collision_name.queryarr("scattering_processes", scattering_process_names);

    // create a vector of ScatteringProcess objects from each scattering
    // process name
    amrex::Vector<ScatteringProcess> scattering_processes;
    for (const auto& scattering_process : scattering_process_names) {
        const std::string kw_cross_section = scattering_process + "_cross_section";
        std::string cross_section_file;
        pp_collision_name.query(kw_cross_section.c_str(), cross_section_file);

        // if the scattering process is excitation or ionization get the
        // energy associated with that process
        amrex::ParticleReal energy = 0._prt;
        if (scattering_process.find("excitation") != std::string::npos ||
            scattering_process.find("ionization") != std::string::npos) {
            const std::string kw_energy = scattering_process + "_energy";
            utils::parser::getWithParser(
                pp_collision_name, kw_energy.c_str(), energy);
        }

        ScatteringProcess process(scattering_process, cross_section_file, energy);

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(process.type() != ScatteringProcessType::INVALID,
                                        "Cannot add an unknown scattering process type");

        scattering_processes.push_back(std::move(process));
    }

    m_process_count = static_cast<int>(scattering_processes.size());

    // Store ScatteringProcess::Executor(s) and cross-section array(s).
    // Note that it is necessary to copy over the cross-section array since
    // ScatteringProcess cannot be copied and the cross-section array
    // cannot be stored in the Executor struct as PODVectors can only
    // contain trivially copyable objects.
#ifdef AMREX_USE_GPU
    amrex::Gpu::HostVector<ScatteringProcess::Executor> h_scattering_processes_exe;
    amrex::Vector< amrex::Gpu::HostVector<amrex::ParticleReal> > h_sigmas;
    for (auto const& p : scattering_processes) {
        h_scattering_processes_exe.push_back(p.executor());
        h_sigmas.push_back(p.getSigmaArray());
    }
    m_scattering_processes_exe.resize(m_process_count);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_scattering_processes_exe.begin(),
                        h_scattering_processes_exe.end(), m_scattering_processes_exe.begin());
    m_sigmas_d.resize(m_process_count);
    for (int ii = 0; ii < m_process_count; ii++) {
        m_sigmas_d[ii].resize(h_sigmas[ii].size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_sigmas[ii].begin(),
                            h_sigmas[ii].end(), m_sigmas[ii].begin());
    }
    amrex::Gpu::streamSynchronize();
#else
    for (auto const& p : scattering_processes) {
        m_scattering_processes_exe.push_back(p.executor());
        m_sigmas.push_back(p.getSigmaArray());
    }
#endif
    // Point the cross-section interpolator to the cross-section data stored
    // in this class.
    for (int ii = 0; ii < m_process_count; ii++) {
        m_scattering_processes_exe[ii].m_sigmas_data = m_sigmas[ii].data();
    }
    m_scattering_processes_data = m_scattering_processes_exe.data();
}
