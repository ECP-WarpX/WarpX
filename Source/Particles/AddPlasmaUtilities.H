/* Copyright 2024 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ADDPLASMAUTILITIES_H_
#define WARPX_ADDPLASMAUTILITIES_H_

#include "Initialization/PlasmaInjector.H"

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

bool find_overlap (const amrex::RealBox& tile_real_box, const amrex::RealBox& part_realbox,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo,
                   amrex::RealBox& overlap_realbox, amrex::Box& overlap_box, amrex::IntVect& shifted);

bool find_overlap_flux (const amrex::RealBox& tile_real_box, const amrex::RealBox& part_realbox,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo,
                        const PlasmaInjector& plasma_injector,
                        amrex::RealBox& overlap_realbox, amrex::Box& overlap_box, amrex::IntVect& shifted);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_scale_fac_volume (const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
                                      const amrex::Long pcount) {
    using namespace amrex::literals;
    amrex::Real scale_fac = 0.0_rt;
    if( pcount != 0) {
#if defined(WARPX_DIM_3D)
        scale_fac = dx[0]*dx[1]*dx[2]/pcount;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        scale_fac = dx[0]*dx[1]/pcount;
#elif defined(WARPX_DIM_1D_Z)
        scale_fac = dx[0]/pcount;
#endif
    }
    return scale_fac;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_area_weights (const amrex::IntVect& iv, const int normal_axis) {
    int r;
#if defined(WARPX_DIM_3D)
    r = iv[0]*iv[1]*iv[2]/iv[normal_axis];
#elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
    r = iv[0]*iv[1];
    if (normal_axis == 0) { r /= iv[0]; }
    else if (normal_axis == 2) { r /= iv[1]; }
#elif defined(WARPX_DIM_1D_Z)
    r = iv[0];
    if (normal_axis == 2) { r /= iv[0]; }
#endif
    return r;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_scale_fac_area (const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
                                    const amrex::Real num_ppc_real, const int flux_normal_axis) {
    using namespace amrex::literals;
    amrex::Real scale_fac = 0._rt;
    // Scale particle weight by the area of the emitting surface, within one cell
#if defined(WARPX_DIM_3D)
    scale_fac = dx[0]*dx[1]*dx[2]/dx[flux_normal_axis]/num_ppc_real;
#elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
    scale_fac = dx[0]*dx[1]/num_ppc_real;
    // When emission is in the r direction, the emitting surface is a cylinder.
    // The factor 2*pi*r is added later below.
    if (flux_normal_axis == 0) { scale_fac /= dx[0]; }
    // When emission is in the z direction, the emitting surface is an annulus
    // The factor 2*pi*r is added later below.
    if (flux_normal_axis == 2) { scale_fac /= dx[1]; }
    // When emission is in the theta direction (flux_normal_axis == 1),
    // the emitting surface is a rectangle, within the plane of the simulation
#elif defined(WARPX_DIM_1D_Z)
    scale_fac = dx[0]/num_ppc_real;
    if (flux_normal_axis == 2) { scale_fac /= dx[0]; }
#endif
    return scale_fac;
}

#endif /*WARPX_ADDPLASMAUTILITIES_H_*/
