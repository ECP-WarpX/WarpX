/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PROTON_BORON_FUSION_CROSS_SECTION_H
#define PROTON_BORON_FUSION_CROSS_SECTION_H

#include <AMReX_REAL.H>


/**
 * \brief Computes the total proton-boron fusion cross section. Needs to be implemented.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::ParticleReal ProtonBoronFusionCrossSection (const amrex::ParticleReal& /*E_star*/)
{
    return amrex::ParticleReal(0.);
}


#endif // PROTON_BORON_FUSION_CROSS_SECTION_H
