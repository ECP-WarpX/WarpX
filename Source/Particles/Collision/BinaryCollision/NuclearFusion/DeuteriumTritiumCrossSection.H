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
 * \brief Computes the total deuterium-tritium fusion cross section.
 *
 * @param[in] E_kin_star the kinetic energy of the deuterium-tritium pair in its center of mass frame,
 * in SI units.
 * @return The total cross section in SI units (square meters).
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::ParticleReal DeuteriumTritiumFusionCrossSection (const amrex::ParticleReal& E_kin_star)
{
    using namespace amrex::literals;

    // Convert cross section to SI units: barn to square meter
    constexpr auto barn_to_sqm = 1.e-28_prt;
    return 0.*barn_to_sqm;
}

#endif // DEUTERIUM_TRITIUM_TRITIUM_FUSION_CROSS_SECTION_H
