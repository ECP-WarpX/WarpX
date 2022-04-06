/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDQUADRUPOLE_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDQUADRUPOLE_H_

#include "HardEdged_K.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>

#include <string>

struct HardEdgedQuadrupole
{

    HardEdgedQuadrupole ();

    ~HardEdgedQuadrupole () = default;

    int nelements = 0;

    amrex::Gpu::DeviceVector<amrex::Real> d_zs;
    amrex::Gpu::DeviceVector<amrex::Real> d_ze;
    amrex::Gpu::DeviceVector<amrex::Real> d_dEdx;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void get_field (const int ielement,
                    const amrex::ParticleReal x,
                    const amrex::ParticleReal y,
                    const amrex::ParticleReal z,
                    const amrex::ParticleReal ux,
                    const amrex::ParticleReal uy,
                    const amrex::ParticleReal uz,
                    const amrex::ParticleReal dt,
                    amrex::ParticleReal& Ex,
                    amrex::ParticleReal& Ey)
    {

        amrex::Real const frac = hard_edged_fraction(z, ux, uy, uz, dt, d_zs[ielement], d_ze[ielement]);

        amrex::Real dEdx = frac*d_dEdx[ielement];

        Ex = +x*dEdx;
        Ey = -y*dEdx;

    }

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDQUADRUPOLE_H_
