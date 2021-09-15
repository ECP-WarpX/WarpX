/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PROTON_BORON_FUSION_CROSS_SECTION_H
#define PROTON_BORON_FUSION_CROSS_SECTION_H

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>

#include <cmath>

/**
 * \brief Computes the total proton-boron fusion cross section using the analytical fits given in
 * W.M. Nevins and R. Swain, Nuclear Fusion, 40, 865 (2000). The result is returned in SI units
 * (square meters).
 * For the record, note that there is a typo in equation (1) of this paper: the total cross section
 * should read S(E)/E*exp(-sqrt(E_G/G)) instead of S(E)/E*exp(sqrt(E_G/G)) (minus sign in the
 * exponential).
 *
 * @param[in] E_kin_star the kinetic energy of the proton-boron pair in its center of mass frame,
 * in SI units.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::ParticleReal ProtonBoronFusionCrossSection (const amrex::ParticleReal& E_kin_star)
{
    // If kinetic energy is 0, return a 0 corss section and avoid later division by 0.
    if (E_kin_star == amrex::ParticleReal(0.)) {return amrex::ParticleReal(0.);}

    // Fits use energy in keV and MeV
    constexpr amrex::ParticleReal joule_to_kev = amrex::ParticleReal(1.e-3)/PhysConst::q_e;
    constexpr amrex::ParticleReal joule_to_mev = amrex::ParticleReal(1.e-6)/PhysConst::q_e;
    const amrex::ParticleReal E_kev = E_kin_star*joule_to_kev;
    const amrex::ParticleReal E_mev = E_kin_star*joule_to_mev;

    // Compute Gamow factor, in MeV
    constexpr auto one_pr = amrex::ParticleReal(1.);
    constexpr auto Z_boron = amrex::ParticleReal(5.);
    constexpr amrex::ParticleReal m_boron = amrex::ParticleReal(10.7319) * PhysConst::m_p;
    constexpr amrex::ParticleReal m_reduced = m_boron / (one_pr + m_boron/PhysConst::m_p);
    constexpr amrex::ParticleReal gamow_factor = m_reduced / amrex::ParticleReal(2.) *
                        (PhysConst::q_e*PhysConst::q_e * Z_boron /
                           (amrex::ParticleReal(2.)*PhysConst::ep0*PhysConst::hbar)) *
                        (PhysConst::q_e*PhysConst::q_e * Z_boron /
                           (amrex::ParticleReal(2.)*PhysConst::ep0*PhysConst::hbar)) *
                        joule_to_mev;

    // Compute astrophysical factor, in MeV barn, using the fits
    constexpr auto E_lim1 = amrex::ParticleReal(400.); // Limits between the different fit regions
    constexpr auto E_lim2 = amrex::ParticleReal(642.);
    amrex::ParticleReal astrophysical_factor;
    if (E_kev < E_lim1)
    {
        constexpr auto C0 = amrex::ParticleReal(197.);
        constexpr auto C1 = amrex::ParticleReal(0.24);
        constexpr auto C2 = amrex::ParticleReal(2.31e-4);
        constexpr auto AL = amrex::ParticleReal(1.82e4);
        constexpr auto EL = amrex::ParticleReal(148.);
        constexpr auto dEL_sq = amrex::ParticleReal(2.35*2.35);
        astrophysical_factor = C0 + C1*E_kev + C2*E_kev*E_kev +
                                               AL/((E_kev - EL)*(E_kev - EL) + dEL_sq);
    }
    else if (E_kev < E_lim2)
    {
        constexpr auto D0 = amrex::ParticleReal(330.);
        constexpr auto D1 = amrex::ParticleReal(66.1);
        constexpr auto D2 = amrex::ParticleReal(-20.3);
        constexpr auto D5 = amrex::ParticleReal(-1.58);
        const amrex::ParticleReal E_norm = (E_kev-amrex::ParticleReal(400.))
                                                            /amrex::ParticleReal(100.);
        astrophysical_factor = D0 + D1*E_norm + D2*E_norm*E_norm + D5*std::pow(E_norm,5);
    }
    else
    {
        constexpr auto A0 = amrex::ParticleReal(2.57e6);
        constexpr auto A1 = amrex::ParticleReal(5.67e5);
        constexpr auto A2 = amrex::ParticleReal(1.34e5);
        constexpr auto A3 = amrex::ParticleReal(5.68e5);
        constexpr auto E0 = amrex::ParticleReal(581.3);
        constexpr auto E1 = amrex::ParticleReal(1083.);
        constexpr auto E2 = amrex::ParticleReal(2405.);
        constexpr auto E3 = amrex::ParticleReal(3344.);
        constexpr auto dE0_sq = amrex::ParticleReal(85.7*85.7);
        constexpr auto dE1_sq = amrex::ParticleReal(234.*234.);
        constexpr auto dE2_sq = amrex::ParticleReal(138.*138.);
        constexpr auto dE3_sq = amrex::ParticleReal(309.*309.);
        constexpr auto B = amrex::ParticleReal(4.38);
        astrophysical_factor =  A0 / ((E_kev-E0)*(E_kev-E0) + dE0_sq) +
                                A1 / ((E_kev-E1)*(E_kev-E1) + dE1_sq) +
                                A2 / ((E_kev-E2)*(E_kev-E2) + dE2_sq) +
                                A3 / ((E_kev-E3)*(E_kev-E3) + dE3_sq) + B;
    }

    // Compute cross section, in barn
    const amrex::ParticleReal cross_section_b = astrophysical_factor/E_mev*
                                                std::exp(-std::sqrt(gamow_factor/E_mev));

    // Convert cross section to SI units: barn to square meter
    constexpr auto barn_to_sqm = amrex::ParticleReal(1.e-28);
    return cross_section_b*barn_to_sqm;
}


#endif // PROTON_BORON_FUSION_CROSS_SECTION_H
