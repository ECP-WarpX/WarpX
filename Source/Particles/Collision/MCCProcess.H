/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_MCCPROCESS_H_
#define WARPX_PARTICLES_COLLISION_MCCPROCESS_H_

#include <AMReX_Math.H>
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

    MCCProcess (MCCProcess const&) = delete;
    MCCProcess& operator= (MCCProcess const&) = delete;

    MCCProcess (MCCProcess &&) = default;
    MCCProcess& operator= (MCCProcess &&) = default;

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
                               amrex::Vector<amrex::Real>& energies,
                               amrex::Gpu::HostVector<amrex::Real>& sigmas
                               );

    static
    void sanityCheckEnergyGrid (
                                const amrex::Vector<amrex::Real>& energies,
                                amrex::Real dE
                                );

    struct Executor {
        /** Get the collision cross-section using a simple linear interpolator. If
         * the energy value is lower (higher) than the given energy range the
         * first (last) cross-section value is used.
         *
         * @param E_coll collision energy in eV
         *
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::Real getCrossSection (amrex::Real E_coll) const
        {
            if (E_coll < m_energy_lo) {
                return m_sigma_lo;
            } else if (E_coll > m_energy_hi) {
                return m_sigma_hi;
            } else {
                using amrex::Math::floor;
                using amrex::Math::ceil;
                // calculate index of bounding energy pairs
                amrex::Real temp = (E_coll - m_energy_lo) / m_dE;
                int idx_1 = static_cast<int>(floor(temp));
                int idx_2 = static_cast<int>(ceil(temp));

                // linearly interpolate to the given energy value
                temp -= idx_1;
                return m_sigmas_data[idx_1] + (m_sigmas_data[idx_2] - m_sigmas_data[idx_1]) * temp;
            }
        }

        amrex::Real* m_sigmas_data = nullptr;
        amrex::Real m_energy_lo, m_energy_hi, m_sigma_lo, m_sigma_hi, m_dE;
        amrex::Real m_energy_penalty;
        MCCProcessType m_type;
    };

    Executor const& executor () const {
#ifdef AMREX_USE_GPU
        return m_exe_d;
#else
        return m_exe_h;
#endif
    }

    amrex::Real getCrossSection (amrex::Real E_coll) const
    {
        return m_exe_h.getCrossSection(E_coll);
    }

    amrex::Real getEnergyPenalty () const { return m_exe_h.m_energy_penalty; }

    MCCProcessType type () const { return m_exe_h.m_type; }

private:

    static
    MCCProcessType parseProcessType(const std::string& process);

    void init (const std::string& scattering_process, const amrex::Real energy);

    amrex::Vector<amrex::Real> m_energies;

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::Real> m_sigmas_d;
    Executor m_exe_d;
#endif
    amrex::Gpu::HostVector<amrex::Real> m_sigmas_h;
    Executor m_exe_h;

    int m_grid_size;
};

#endif // WARPX_PARTICLES_COLLISION_MCCPROCESS_H_
