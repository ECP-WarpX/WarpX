/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef LINEAR_BREIT_WHEELER_CROSS_SECTION_H
#define LINEAR_BREIT_WHEELER_CROSS_SECTION_H

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>

#include <cmath>

/**
 * \brief Computes the total proton-boron fusion cross section. When E_kin_star < 3.5 MeV, we use
 * the analytical fits given in W.M. Nevins and R. Swain, Nuclear Fusion, 40, 865 (2000). When
 * E_kin_star > 3.5 MeV, we use a simple power law fit of the data presented in Buck et al.,
 * Nuclear Physics A, 398(2), 189-202 (1983). Both fits return the same value for
 * E_kin_star = 3.5 MeV.
 *
 * @param[in] E_kin_star the kinetic energy of the proton-boron pair in its center of mass frame,
 * in SI units.
 * @return The total cross section in SI units (square meters).
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::ParticleReal LinearBreitWheelerCrossSection (const amrex::ParticleReal& E_star)
{
    using namespace amrex::literals;

    constexpr amrex::ParticleReal c_sq = PhysConst::c * PhysConst::c;
    constexpr amrex::ParticleReal me = PhysConst::m_e;
    // If kinetic energy is 0, return a 0 cross section and avoid later division by 0.
    if (E_star < 2._prt*me*c_sq) {return 0._prt;}

    constexpr auto one_half_pr = amrex::ParticleReal(1./2.);
    constexpr amrex::ParticleReal pi = MathConst::pi;
    constexpr amrex::ParticleReal re2 = PhysConst::r_e * PhysConst::r_e;
    auto constexpr pow2 = [](double const x) { return x*x; };
    auto constexpr pow4 = [](double const x) { return x*x*x*x; };

    const amrex::ParticleReal s = E_star/(me * c_sq);
    const amrex::ParticleReal beta = std::sqrt(1._prt-1._prt/s);

    const amrex::ParticleReal term1 = (3._prt-pow4(beta))*std::log((1._prt+beta)/(1._prt-beta));
    const amrex::ParticleReal term2 = 2._prt*beta*(pow2(beta)-2._prt);
    const amrex::ParticleReal cross_section = one_half_pr*pi*re2*(1._prt-pow2(beta))*(term1+term2);
    
    //amrex::Print() <<  "PRINTING CROSS SECTION" << cross_section << "\n";  

    return cross_section;
}

#endif // LINEAR_BREIT_WHEELER_CROSS_SECTION_H
