/* Copyright 2019 Maxence Thevenet, Remi Lehe, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SUM_GUARD_CELLS_H_
#define WARPX_SUM_GUARD_CELLS_H_

#include <AMReX_MultiFab.H>

/** \brief Sum the values of `mf`, where the different boxes overlap
 * (i.e. in the guard cells)
 *
 * This is typically called for the sources of the Maxwell equations (J/rho)
 * after deposition from the macroparticles.
 *
 *  - When WarpX is used with a finite-difference scheme: this only
 *    updates the *valid* cells of `mf`
 *  - When WarpX is used with a spectral scheme (PSATD): this
 *    updates both the *valid* cells and *guard* cells. (This is because a
 *    spectral solver requires the value of the sources over a large stencil.)
 */
void
WarpXSumGuardCells(amrex::MultiFab& mf, const amrex::Periodicity& period,
                   const amrex::IntVect& src_ngrow,
                   const int icomp=0, const int ncomp=1);

/** \brief Sum the values of `src` where the different boxes overlap
 * (i.e. in the guard cells) and copy them into `dst`
 *
 * This is typically called for the sources of the Maxwell equations (J/rho)
 * after deposition from the macroparticles + filtering.
 *
 *  - When WarpX is used with a finite-difference scheme: this only
 *    updates the *valid* cells of `dst`
 *  - When WarpX is used with a spectral scheme (PSATD): this
 *    updates both the *valid* cells and *guard* cells. (This is because a
 *    spectral solver requires the value of the sources over a large stencil.)
 *
 * Note: `i_comp` is the component where the results will be stored in `dst`;
 *       The component from which we copy in `src` is always 0.
 */
void
WarpXSumGuardCells(amrex::MultiFab& dst, amrex::MultiFab& src,
                   const amrex::Periodicity& period,
                   const amrex::IntVect& src_ngrow,
                   const int icomp=0, const int ncomp=1);

#endif // WARPX_SUM_GUARD_CELLS_H_
