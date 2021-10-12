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
    using namespace amrex::literals;

    // If kinetic energy is 0, return a 0 cross section and avoid later division by 0.
    if (E_kin_star == 0._prt) {return 0._prt;}

    // Fits use energy in keV and MeV
    constexpr amrex::ParticleReal joule_to_kev = 1.e-3_prt/PhysConst::q_e;
    constexpr amrex::ParticleReal joule_to_mev = 1.e-6_prt/PhysConst::q_e;
    const amrex::ParticleReal E_kev = E_kin_star*joule_to_kev;
    const amrex::ParticleReal E_mev = E_kin_star*joule_to_mev;

    // Compute Gamow factor, in MeV
    constexpr auto one_pr = 1._prt;
    constexpr auto Z_boron = 5._prt;
    constexpr amrex::ParticleReal m_boron = 10.7319_prt * PhysConst::m_p;
    constexpr amrex::ParticleReal m_reduced = m_boron / (one_pr + m_boron/PhysConst::m_p);
    constexpr amrex::ParticleReal gamow_factor = m_reduced / 2._prt *
                        (PhysConst::q_e*PhysConst::q_e * Z_boron /
                           (2._prt*PhysConst::ep0*PhysConst::hbar)) *
                        (PhysConst::q_e*PhysConst::q_e * Z_boron /
                           (2._prt*PhysConst::ep0*PhysConst::hbar)) *
                        joule_to_mev;

    // Compute astrophysical factor, in MeV barn, using the fits
    constexpr auto E_lim1 = 400._prt; // Limits between the different fit regions
    constexpr auto E_lim2 = 642._prt;
    amrex::ParticleReal astrophysical_factor;
    if (E_kev < E_lim1)
    {
        constexpr auto C0 = 197._prt;
        constexpr auto C1 = 0.24_prt;
        constexpr auto C2 = 2.31e-4_prt;
        constexpr auto AL = 1.82e4_prt;
        constexpr auto EL = 148._prt;
        constexpr auto dEL_sq = 2.35_prt*2.35_prt;
        astrophysical_factor = C0 + C1*E_kev + C2*E_kev*E_kev +
                                               AL/((E_kev - EL)*(E_kev - EL) + dEL_sq);
    }
    else if (E_kev < E_lim2)
    {
        constexpr auto D0 = 330._prt;
        constexpr auto D1 = 66.1_prt;
        constexpr auto D2 = -20.3_prt;
        constexpr auto D5 = -1.58_prt;
        const amrex::ParticleReal E_norm = (E_kev-400._prt) * 1.e-2_prt;
        astrophysical_factor = D0 + D1*E_norm + D2*E_norm*E_norm + D5*std::pow(E_norm,5);
    }
    else
    {
        constexpr auto A0 = 2.57e6_prt;
        constexpr auto A1 = 5.67e5_prt;
        constexpr auto A2 = 1.34e5_prt;
        constexpr auto A3 = 5.68e5_prt;
        constexpr auto E0 = 581.3_prt;
        constexpr auto E1 = 1083._prt;
        constexpr auto E2 = 2405._prt;
        constexpr auto E3 = 3344._prt;
        constexpr auto dE0_sq = 85.7_prt*85.7_prt;
        constexpr auto dE1_sq = 234._prt*234._prt;
        constexpr auto dE2_sq = 138._prt*138._prt;
        constexpr auto dE3_sq = 309._prt*309._prt;
        constexpr auto B = 4.38_prt;
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
