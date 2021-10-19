/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PARTICLES_KINETIC_ENERGY_H_
#define PARTICLES_KINETIC_ENERGY_H_

#include "Utils/WarpXConst.H"

#include "AMReX_Extension.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_REAL.H"

#include <cmath>

namespace Algorithms{

    // This marks the gamma threshold to switch between the full relativistic expression
    // for particle kinetic energy and a Taylor expansion.
    static constexpr auto gamma_relativistic_threshold =
        static_cast<amrex::Real>(1.005);

    /**
     * \brief Computes the kinetic energy of a particle. Below a threshold for the
     * Lorentz factor (gamma_relativistic_threshold) it uses a Taylor expansion instead of
     * the full relativistic expression. This method should not be used with photons.
     *
     * @param[in] ux x component of the particle momentum (code units)
     * @param[in] uy y component of the particle momentum (code units)
     * @param[in] uz z component of the particle momentum (code units)
     * @param[in] mass mass of the particle (in S.I. units)
     *
     * @return the kinetic energy of the particle (in S.I. units)
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real KineticEnergy(
        const amrex::Real ux, const amrex::Real uy, const amrex::Real uz,
        const amrex::Real mass)
    {
        using namespace amrex;

        constexpr auto c2 = PhysConst::c * PhysConst::c;
        constexpr auto inv_c2 = 1.0_rt/c2;

        const auto u2 = (ux*ux + uy*uy + uz*uz)*inv_c2;
        const auto gamma = std::sqrt(1.0_rt + u2);

        const auto kk = (gamma > gamma_relativistic_threshold)?
            (gamma-1.0_rt):
            (u2*0.5_rt - u2*u2*(1.0_rt/8.0_rt) + u2*u2*u2*(1.0_rt/16.0_rt)-
            u2*u2*u2*u2*(5.0_rt/128.0_rt) + (7.0_rt/256_rt)*u2*u2*u2*u2*u2); //Taylor expansion

        return kk*mass*c2;
    }

    /**
     * \brief Computes the kinetic energy of a photon.
     *
     * @param[in] ux x component of the particle momentum (code units)
     * @param[in] uy y component of the particle momentum (code units)
     * @param[in] uz z component of the particle momentum (code units)
     *
     * @return the kinetic energy of the photon (in S.I. units)
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real KineticEnergyPhotons(
        const amrex::Real ux, const amrex::Real uy, const amrex::Real uz)
    {
        // Photons have zero mass, but ux, uy and uz are calculated assuming a mass equal to the
        // electron mass. Hence, photons need a special treatment to calculate the total energy.
        constexpr auto me_c = PhysConst::m_e * PhysConst::c;

        return me_c * std::sqrt(ux*ux + uy*uy + uz*uz);
    }

}

#endif // PARTICLES_ALGORITHMS_H_
