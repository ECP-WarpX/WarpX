/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDPLASMALENS_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDPLASMALENS_H_

#include "LatticeElementBase.H"
#include "HardEdged_K.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>

#include <string>

struct HardEdgedPlasmaLensDevice;

struct HardEdgedPlasmaLens
    : LatticeElementBase
{

    HardEdgedPlasmaLens ();

    ~HardEdgedPlasmaLens () = default;

    amrex::Gpu::DeviceVector<amrex::Real> d_dEdx;
    amrex::Gpu::DeviceVector<amrex::Real> d_dBdx;

    HardEdgedPlasmaLensDevice GetDeviceInstance () const;


};

struct HardEdgedPlasmaLensDevice
{

    void InitHardEdgedPlasmaLensDevice (HardEdgedPlasmaLens const& h_plasmalens);

    int nelements = 0;

    const amrex::Real* AMREX_RESTRICT d_zs_arr;
    const amrex::Real* AMREX_RESTRICT d_ze_arr;

    const amrex::Real* AMREX_RESTRICT d_dEdx_arr;
    const amrex::Real* AMREX_RESTRICT d_dBdx_arr;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void get_field (const int ielement,
                    const amrex::ParticleReal x,
                    const amrex::ParticleReal y,
                    const amrex::ParticleReal z,
                    const amrex::ParticleReal zpvdt,
                    amrex::ParticleReal& Ex,
                    amrex::ParticleReal& Ey,
                    amrex::ParticleReal& Bx,
                    amrex::ParticleReal& By) const
    {

        amrex::Real const frac = hard_edged_fraction(z, zpvdt, d_zs_arr[ielement], d_ze_arr[ielement]);

        amrex::Real dEdx = frac*d_dEdx_arr[ielement];
        amrex::Real dBdx = frac*d_dBdx_arr[ielement];

        Ex = +x*dEdx;
        Ey = +y*dEdx;
        Bx = +y*dBdx;
        By = -x*dBdx;

    }

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDPLASMALENS_H_
