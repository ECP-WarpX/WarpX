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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real hard_edged_fraction(const amrex::Real z,
                                const amrex::Real ux,
                                const amrex::Real uy,
                                const amrex::Real uz,
                                const amrex::Real dt,
                                const amrex::Real zs,
                                const amrex::Real ze)
{

    using namespace amrex::literals;

    constexpr amrex::Real inv_c2 = 1._rt/(PhysConst::c*PhysConst::c);
    const amrex::ParticleReal gamma = std::sqrt(1._rt + (ux*ux + uy*uy + uz*uz)*inv_c2);
    const amrex::ParticleReal vz = uz/gamma;

    amrex::ParticleReal const z1 = z;
    amrex::ParticleReal const z2 = z + vz*dt;

    // This allows vz to be positive or negative
    amrex::ParticleReal const zl = std::min(z1, z2);
    amrex::ParticleReal const zr = std::max(z1, z2);

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
