/* Copyright 2019 Remi Lehe, Revathi Jambunathan, Revathi Jambunathan
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PML_KERNELS_H_
#define WARPX_PML_KERNELS_H_

#include <AMReX.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ex (int i, int j, int k, Array4<Real> const& Ex,
                        const Real* const sigma_fac_y,
                        const Real* const sigma_fac_z,
                        const Real* const sigma_star_fac_x,
                        int xlo,int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
    Ex(i,j,k,0) *= sigma_fac_y[j-ylo];
    Ex(i,j,k,1) *= sigma_fac_z[k-zlo];
#else
    Ex(i,j,k,1) *= sigma_fac_z[j-zlo];
    amrex::ignore_unused(sigma_fac_y, ylo);
#endif
    Ex(i,j,k,2) *= sigma_star_fac_x[i-xlo];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ey (int i, int j, int k, Array4<Real> const& Ey,
                        const Real* const sigma_fac_z,
                        const Real* const sigma_fac_x,
                        const Real* const sigma_star_fac_y,
                        int xlo,int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
    Ey(i,j,k,0) *= sigma_fac_z[k-zlo];
    Ey(i,j,k,2) *= sigma_star_fac_y[j-ylo];
#else
    Ey(i,j,k,0) *= sigma_fac_z[j-zlo];
    amrex::ignore_unused(sigma_star_fac_y, ylo);
#endif
    Ey(i,j,k,1) *= sigma_fac_x[i-xlo];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ez (int i, int j, int k, Array4<Real> const& Ez,
                        const Real* const sigma_fac_x,
                        const Real* const sigma_fac_y,
                        const Real* const sigma_star_fac_z,
                        int xlo,int ylo, int zlo)
{
    Ez(i,j,k,0) *= sigma_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
    Ez(i,j,k,1) *= sigma_fac_y[j-ylo];
    Ez(i,j,k,2) *= sigma_star_fac_z[k-zlo];
#else
    Ez(i,j,k,2) *= sigma_star_fac_z[j-zlo];
    amrex::ignore_unused(sigma_fac_y, ylo);
#endif
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_bx (int i, int j, int k, Array4<Real> const& Bx,
                        const Real* const sigma_star_fac_y,
                        const Real* const sigma_star_fac_z,
                        int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
   Bx(i,j,k,0) *= sigma_star_fac_y[j-ylo];
   Bx(i,j,k,1) *= sigma_star_fac_z[k-zlo];
#else
   Bx(i,j,k,1) *= sigma_star_fac_z[j-zlo];
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
   By(i,j,k,0) *= sigma_star_fac_z[k-zlo];
#else
   By(i,j,k,0) *= sigma_star_fac_z[j-zlo];
#endif
   By(i,j,k,1) *= sigma_star_fac_x[i-xlo];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_bz (int i, int j, int k, Array4<Real> const& Bz,
                        const Real* const sigma_star_fac_x,
                        const Real* const sigma_star_fac_y,
                        int xlo, int ylo)
{
   Bz(i,j,k,0) *= sigma_star_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
   Bz(i,j,k,1) *= sigma_star_fac_y[j-ylo];
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
   F_fab(i,j,k,0) *= sigma_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
   F_fab(i,j,k,1) *= sigma_fac_y[j-ylo];
   F_fab(i,j,k,2) *= sigma_fac_z[k-zlo];
#else
   F_fab(i,j,k,2) *= sigma_fac_z[j-zlo];
   amrex::ignore_unused(sigma_fac_y, ylo);
#endif
}

#endif
