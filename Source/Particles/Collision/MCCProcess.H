/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_MCCPROCESS_H_
#define WARPX_PARTICLES_COLLISION_MCCPROCESS_H_

#include <AMReX_Vector.H>
#include <AMReX_RandomEngine.H>
#include <AMReX_GpuContainers.H>

enum class MCCProcessType {
    INVALID,
    ELASTIC,
    BACK,
    CHARGE_EXCHANGE,
    EXCITATION,
    IONIZATION,
};

class MCCProcess
{
public:
    template <typename T>
    using VectorType = amrex::Gpu::ManagedVector<T>;

    MCCProcess (
                const std::string& scattering_process,
                const std::string& cross_section_file,
                const amrex::Real energy
                );

    template <typename InputVector>
    MCCProcess (
                const std::string& scattering_process,
                const InputVector&& energies,
                const InputVector&& sigmas,
                const amrex::Real energy
                );

    void* operator new  ( std::size_t count );
    void operator delete (void* ptr);

    /** Get the collision cross-section using a simple linear interpolator. If
     * the energy value is lower (higher) than the given energy range the
     * first (last) cross-section value is used.
     *
     * @param E_coll collision energy in eV
     *
     */
    AMREX_GPU_HOST_DEVICE
    amrex::Real getCrossSection (amrex::Real E_coll) const
    {
        if (E_coll < m_energy_lo) {
            return m_sigma_lo;
        } else if (E_coll > m_energy_hi) {
            return m_sigma_hi;
        } else {
            using std::floor;
            using std::ceil;
            // calculate index of bounding energy pairs
            amrex::Real temp = (E_coll - m_energy_lo) / m_dE;
            int idx_1 = floor(temp);
            int idx_2 = ceil(temp);

            // linearly interpolate to the given energy value
            temp -= idx_1;
            return (
                    m_sigmas_data[idx_1] + (m_sigmas_data[idx_2] - m_sigmas_data[idx_1]) * temp
                    );
        }
    }

    /** Read the given cross-section data file to memory.
     *
     * @param cross_section_file the path to the file containing the cross-
     *        section data
     * @param energies vector storing energy values in eV
     * @param sigmas vector storing cross-section values
     *
     */
    static
    void readCrossSectionFile (
                               const std::string cross_section_file,
                               VectorType<amrex::Real>& energies,
                               VectorType<amrex::Real>& sigmas
                               );

    static
    void sanityCheckEnergyGrid (
                                const VectorType<amrex::Real>& energies,
                                amrex::Real dE
                                );

    MCCProcessType m_type;
    amrex::Real m_energy_penalty;

private:

    static
    MCCProcessType parseProcessType(const std::string& process);

    void init();

    VectorType<amrex::Real> m_energies;
    amrex::Real* m_energies_data;
    VectorType<amrex::Real> m_sigmas;
    amrex::Real* m_sigmas_data;

    int m_grid_size;
    amrex::Real m_energy_lo, m_energy_hi, m_sigma_lo, m_sigma_hi, m_dE;
};

#endif // WARPX_PARTICLES_COLLISION_MCCPROCESS_H_
