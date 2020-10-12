/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleNumber.H"
#include "WarpX.H"

// constructor
ParticleNumber::ParticleNumber (std::string rd_name)
: ReducedDiags{rd_name}
{
    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    const auto & mypc = warpx.GetPartContainer();

    // get number of species (int)
    const auto nSpecies = mypc.nSpecies();

    // resize data array to nSpecies+1
    // (number of particles of each species + total number of particles)
    m_data.resize(nSpecies+1, amrex::Real(0.));

    // get species names (std::vector<std::string>)
    const auto species_names = mypc.GetSpeciesNames();

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app};
            // write header row
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            ofs << m_sep;
            ofs << "[3]total()";
            constexpr int shift_first_species = 4; // Column number of first species in output file
            for (int i = 0; i < nSpecies; ++i)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_first_species+i) + "]";
                ofs << species_names[i]+"()";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes total number of macroparticles
void ParticleNumber::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get MultiParticleContainer class object
    const auto & mypc = WarpX::GetInstance().GetPartContainer();

    // get number of species (int)
    const auto nSpecies = mypc.nSpecies();

    // Index of total number of particles (all species) in m_data
    constexpr int idx_total = 0;
    // Index of first species in m_data
    constexpr int idx_first_species = 1;

    // Initialize total number of particles (all species) to 0
    m_data[idx_total] = amrex::Real(0.);

    // loop over species
    for (int i_s = 0; i_s < nSpecies; ++i_s)
    {
        // get WarpXParticleContainer class object
        const auto & myspc = mypc.GetParticleContainer(i_s);

        // Save result for this species
        m_data[idx_first_species + i_s] = myspc.TotalNumberOfParticles();
        // Increase total number of particles (all species)
        m_data[idx_total] += m_data[idx_first_species + i_s];
    }
    // end loop over species

    /* m_data now contains up-to-date values for:
     *  [total number of macroparticles (all species),
     *   total number of macroparticles (species 1),
     *   ...,
     *   total number of macroparticles (species n)] */

}
// end void ParticleNumber::ComputeDiags
