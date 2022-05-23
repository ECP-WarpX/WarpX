/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGED_K_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGED_K_H_

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * \brief Calculate the residence correction, the fraction of the time step the particle
 * spends inside of the hard edge lattice element
 *
 * @param[in] z the current location of the particle
 * @param[in] zpvdt the estimated future location of the particle, z + v*dt
 * @param[in] zs the start of the lattice element
 * @param[in] ze the end of the lattice element
 * @param[out] the fraction is returned
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real hard_edged_fraction(const amrex::Real z,
                                const amrex::Real zpvdt,
                                const amrex::Real zs,
                                const amrex::Real ze)
{

    using namespace amrex::literals;

    // This allows vz to be positive or negative
    amrex::ParticleReal const zl = std::min(z, zpvdt);
    amrex::ParticleReal const zr = std::max(z, zpvdt);

    // Calculate the residence correction
    // frac will be 1 if the step is completely inside the lens, between 0 and 1
    // when entering or leaving the lens, and otherwise 0.
    amrex::Real const fl = ( (zl >= zs && zl < ze) ? 1. : 0.);
    amrex::Real const fr = ( (zr >= zs && zr < ze) ? 1. : 0.);
    amrex::Real frac = fl;
    if (fl > fr) frac = (ze - zl)/(zr - zl);
    if (fr > fl) frac = (zr - zs)/(zr - zl);

    return frac;
}

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGED_K_H_
