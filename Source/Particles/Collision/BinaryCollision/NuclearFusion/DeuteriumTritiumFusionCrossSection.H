/* Copyright 2022 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef DEUTERIUM_TRITIUM_FUSION_CROSS_SECTION_H
#define DEUTERIUM_TRITIUM_FUSION_CROSS_SECTION_H

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>

#include <cmath>

/**
 * \brief Computes the total deuterium-tritium fusion cross section, using
 * the analytical fits given in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611
 *
 * @param[in] E_kin_star the kinetic energy of the deuterium-tritium pair in its center of mass frame, in SI units.
 * @return The total cross section in SI units (square meters).
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::ParticleReal DeuteriumTritiumFusionCrossSection (const amrex::ParticleReal& E_kin_star)
{
    using namespace amrex::literals;

    constexpr amrex::ParticleReal joule_to_keV = 1.e-3_prt/PhysConst::q_e;
    const amrex::ParticleReal E_keV = E_kin_star*joule_to_keV;

    // If kinetic energy is 0, return a 0 cross section and avoid later division by 0.
    if (E_keV == 0._prt) {return 0._prt;}

    // Compute the Gamow constant B_G (in keV^{1/2})
    // (See Eq. 3 in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611)
    constexpr amrex::ParticleReal m_D = 2.01410177812 * PhysConst::m_u;
    constexpr amrex::ParticleReal m_T = 3.0160492779 * PhysConst::m_u;
    constexpr amrex::ParticleReal m_reduced = m_D / (1._prt + m_D/m_T);
    amrex::ParticleReal B_G = MathConst::pi * PhysConst::alpha * 1._prt * 1._prt
         * std::sqrt( 2._prt*m_reduced*PhysConst::c*PhysConst::c * joule_to_keV );

    // Compute astrophysical_factor
    // (See Eq. 9 and Table IV in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611)
    constexpr amrex::ParticleReal A1 = 6.927e4;
    constexpr amrex::ParticleReal A2 = 7.454e8;
    constexpr amrex::ParticleReal A3 = 2.050e6;
    constexpr amrex::ParticleReal A4 = 5.2002e4;
    constexpr amrex::ParticleReal B1 = 6.38e1;
    constexpr amrex::ParticleReal B2 = -9.95e-1;
    constexpr amrex::ParticleReal B3 = 6.981e-5;
    constexpr amrex::ParticleReal B4 = 1.728e-4;

    amrex::ParticleReal astrophysical_factor =
        (A1 + E_keV*(A2 + E_keV*(A3 + E_keV*A4))) /
        (1_prt + E_keV*(B1 + E_keV*(B2 + E_keV*(B3 + E_keV*B4))));

    // Compute cross-section in SI units
    // (See Eq. 8 in H.-S. Bosch and G.M. Hale 1992 Nucl. Fusion 32 611)
    constexpr amrex::ParticleReal millibarn_to_sqm = 1.e-31_prt;
    return millibarn_to_sqm * astrophysical_factor/E_keV * std::exp(-B_G/std::sqrt(E_keV));
}

#endif // DEUTERIUM_TRITIUM_TRITIUM_FUSION_CROSS_SECTION_H
