/* Copyright 2021 Roelof Groenewald
 *
 * This file is part of WarpX.
 *
 * License: ????
 */
#ifndef WARPX_PARTICLES_COLLISION_MCCPROCESS_H_
#define WARPX_PARTICLES_COLLISION_MCCPROCESS_H_

#include <AMReX_Vector.H>

class MCCProcess
{
public:

    MCCProcess (
        const std::string scattering_process,
        const std::string cross_section_file,
        const amrex::Real energy
    );

    /** Get the collision cross-section using a simple linear interpolator. If
     * the energy value is lower (higher) than the given energy range the
     * first (last) cross-section value is used.
     *
     * @param E_coll collision energy in eV
     *
     */
    amrex::Real getCrossSection (amrex::Real E_coll) const;

    /** Read the given cross-section data file to memory.
     *
     * @param cross_section_file the path to the file containing the cross-
     *        section data
     * @param energies vector storing energy values in eV
     * @param sigmas vector storing cross-section values
     *
     */

    void readCrossSectionFile (
        const std::string cross_section_file,
        amrex::Vector<amrex::Real>& energies,
        amrex::Vector<amrex::Real>& sigmas
    );

    void sanityCheckEnergyGrid ( );

    std::string name;
    amrex::Real energy_penalty;

private:

    amrex::Vector<amrex::Real> m_energies;
    amrex::Vector<amrex::Real> m_sigmas;

    int m_grid_size;
    amrex::Real m_energy_lo, m_energy_hi, m_sigma_lo, m_sigma_hi, m_dE;
};

#endif // WARPX_PARTICLES_COLLISION_MCCPROCESS_H_
