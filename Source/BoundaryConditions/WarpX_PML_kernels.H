/* Copyright 2019 Remi Lehe, Revathi Jambunathan, Revathi Jambunathan
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PML_KERNELS_H_
#define WARPX_PML_KERNELS_H_

#include "BoundaryConditions/PMLComponent.H"
#include <AMReX.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ex (int i, int j, int k, Array4<Real> const& Ex,
                        amrex::GpuArray<int, AMREX_SPACEDIM> const& Ex_stag,
                        const Real* const sigma_fac_x,
                        const Real* const sigma_fac_y,
                        const Real* const sigma_fac_z,
                        const Real* const sigma_star_fac_x,
                        const Real* const sigma_star_fac_y,
                        const Real* const sigma_star_fac_z,
                        int xlo, int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 2)

    amrex::ignore_unused(sigma_fac_y, ylo);

    // sx = 0 means that Ex is staggered in x, while sx = 1 means that Ex is nodal in x (same for z)
    const int sx = Ex_stag[0];
    const int sz = Ex_stag[1];

    // Exx
    if (sx == 0) {
        Ex(i,j,k,PMLComp::xx) *= sigma_star_fac_x[i-xlo];
    } else {
        Ex(i,j,k,PMLComp::xx) *= sigma_fac_x[i-xlo];
    }

    // Exz
    if (sz == 0) {
        Ex(i,j,k,PMLComp::xz) *= sigma_star_fac_z[j-zlo];
    } else {
        Ex(i,j,k,PMLComp::xz) *= sigma_fac_z[j-zlo];
    }

#elif (AMREX_SPACEDIM == 3)

    // sx = 0 means that Ex is staggered in x, while sx = 1 means that Ex is nodal in x (same for y, z)
    const int sx = Ex_stag[0];
    const int sy = Ex_stag[1];
    const int sz = Ex_stag[2];

    // Exx
    if (sx == 0) {
        Ex(i,j,k,PMLComp::xx) *= sigma_star_fac_x[i-xlo];
    } else {
        Ex(i,j,k,PMLComp::xx) *= sigma_fac_x[i-xlo];
    }

    // Exy
    if (sy == 0) {
        Ex(i,j,k,PMLComp::xy) *= sigma_star_fac_y[j-ylo];
    } else {
        Ex(i,j,k,PMLComp::xy) *= sigma_fac_y[j-ylo];
    }

    // Exz
    if (sz == 0) {
        Ex(i,j,k,PMLComp::xz) *= sigma_star_fac_z[k-zlo];
    } else {
        Ex(i,j,k,PMLComp::xz) *= sigma_fac_z[k-zlo];
    }

#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ey (int i, int j, int k, Array4<Real> const& Ey,
                        amrex::GpuArray<int, AMREX_SPACEDIM> const& Ey_stag,
                        const Real* const sigma_fac_x,
                        const Real* const sigma_fac_y,
                        const Real* const sigma_fac_z,
                        const Real* const sigma_star_fac_x,
                        const Real* const sigma_star_fac_y,
                        const Real* const sigma_star_fac_z,
                        int xlo, int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 2)

    amrex::ignore_unused(sigma_star_fac_y, ylo);

    // sx = 0 means that Ey is staggered in x, while sx = 1 means that Ey is nodal in x (same for z)
    const int sx = Ey_stag[0];
    const int sz = Ey_stag[1];

    // Eyx
    if (sx == 0) {
        Ey(i,j,k,PMLComp::yx) *= sigma_star_fac_x[i-xlo];
    } else {
        Ey(i,j,k,PMLComp::yx) *= sigma_fac_x[i-xlo];
    }

    // Eyz
    if (sz == 0) {
        Ey(i,j,k,PMLComp::yz) *= sigma_star_fac_z[j-zlo];
    } else {
        Ey(i,j,k,PMLComp::yz) *= sigma_fac_z[j-zlo];
    }

#elif (AMREX_SPACEDIM == 3)

    // sx = 0 means that Ey is staggered in x, while sx = 1 means that Ey is nodal in x (same for y, z)
    const int sx = Ey_stag[0];
    const int sy = Ey_stag[1];
    const int sz = Ey_stag[2];

    // Eyx
    if (sx == 0) {
        Ey(i,j,k,PMLComp::yx) *= sigma_star_fac_x[i-xlo];
    } else {
        Ey(i,j,k,PMLComp::yx) *= sigma_fac_x[i-xlo];
    }

    // Eyy
    if (sy == 0) {
        Ey(i,j,k,PMLComp::yy) *= sigma_star_fac_y[j-ylo];
    } else {
        Ey(i,j,k,PMLComp::yy) *= sigma_fac_y[j-ylo];
    }

    // Eyz
    if (sz == 0) {
        Ey(i,j,k,PMLComp::yz) *= sigma_star_fac_z[k-zlo];
    } else {
        Ey(i,j,k,PMLComp::yz) *= sigma_fac_z[k-zlo];
    }

#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ez (int i, int j, int k, Array4<Real> const& Ez,
                        amrex::GpuArray<int, AMREX_SPACEDIM> const& Ez_stag,
                        const Real* const sigma_fac_x,
                        const Real* const sigma_fac_y,
                        const Real* const sigma_fac_z,
                        const Real* const sigma_star_fac_x,
                        const Real* const sigma_star_fac_y,
                        const Real* const sigma_star_fac_z,
                        int xlo, int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 2)

    amrex::ignore_unused(sigma_fac_y, ylo);

    // sx = 0 means that Ez is staggered in x, while sx = 1 means that Ez is nodal in x (same for z)
    const int sx = Ez_stag[0];
    const int sz = Ez_stag[1];

    // Ezx
    if (sx == 0) {
        Ez(i,j,k,PMLComp::zx) *= sigma_star_fac_x[i-xlo];
    } else {
        Ez(i,j,k,PMLComp::zx) *= sigma_fac_x[i-xlo];
    }

    // Ezz
    if (sz == 0) {
        Ez(i,j,k,PMLComp::zz) *= sigma_star_fac_z[j-zlo];
    } else {
        Ez(i,j,k,PMLComp::zz) *= sigma_fac_z[j-zlo];
    }

#elif (AMREX_SPACEDIM == 3)

    // sx = 0 means that Ez is staggered in x, while sx = 1 means that Ez is nodal in x (same for y, z)
    const int sx = Ez_stag[0];
    const int sy = Ez_stag[1];
    const int sz = Ez_stag[2];

    // Ezx
    if (sx == 0) {
        Ez(i,j,k,PMLComp::zx) *= sigma_star_fac_x[i-xlo];
    } else {
        Ez(i,j,k,PMLComp::zx) *= sigma_fac_x[i-xlo];
    }

    // Ezy
    if (sy == 0) {
        Ez(i,j,k,PMLComp::zy) *= sigma_star_fac_y[j-ylo];
    } else {
        Ez(i,j,k,PMLComp::zy) *= sigma_fac_y[j-ylo];
    }

    // Ezz
    if (sz == 0) {
        Ez(i,j,k,PMLComp::zz) *= sigma_star_fac_z[k-zlo];
    } else {
        Ez(i,j,k,PMLComp::zz) *= sigma_fac_z[k-zlo];
    }

#endif
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_bx (int i, int j, int k, Array4<Real> const& Bx,
                        const Real* const sigma_star_fac_y,
                        const Real* const sigma_star_fac_z,
                        int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
   Bx(i,j,k,PMLComp::xy) *= sigma_star_fac_y[j-ylo];
   Bx(i,j,k,PMLComp::xz) *= sigma_star_fac_z[k-zlo];
#else
   Bx(i,j,k,PMLComp::xz) *= sigma_star_fac_z[j-zlo];
   amrex::ignore_unused(sigma_star_fac_y, ylo);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_by (int i, int j, int k, Array4<Real> const& By,
                        const Real* const sigma_star_fac_z,
                        const Real* const sigma_star_fac_x,
                        int zlo, int xlo)
{
#if (AMREX_SPACEDIM == 3)
   By(i,j,k,PMLComp::yz) *= sigma_star_fac_z[k-zlo];
#else
   By(i,j,k,PMLComp::yz) *= sigma_star_fac_z[j-zlo];
#endif
   By(i,j,k,PMLComp::yx) *= sigma_star_fac_x[i-xlo];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_bz (int i, int j, int k, Array4<Real> const& Bz,
                        const Real* const sigma_star_fac_x,
                        const Real* const sigma_star_fac_y,
                        int xlo, int ylo)
{
   Bz(i,j,k,PMLComp::zx) *= sigma_star_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
   Bz(i,j,k,PMLComp::zy) *= sigma_star_fac_y[j-ylo];
#else
   amrex::ignore_unused(sigma_star_fac_y, ylo);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_F (int i, int j, int k, Array4<Real> const& F_fab,
                      const Real* const sigma_fac_x,
                      const Real* const sigma_fac_y,
                      const Real* const sigma_fac_z,
                      int xlo, int ylo, int zlo)
{
   F_fab(i,j,k,PMLComp::x) *= sigma_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
   F_fab(i,j,k,PMLComp::y) *= sigma_fac_y[j-ylo];
   F_fab(i,j,k,PMLComp::z) *= sigma_fac_z[k-zlo];
#else
   F_fab(i,j,k,PMLComp::z) *= sigma_fac_z[j-zlo];
   amrex::ignore_unused(sigma_fac_y, ylo);
#endif
}

#endif
