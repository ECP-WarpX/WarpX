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

/** \brief Sum overlapping cell values in `mf` independently for each
 * selected component.
 *
 * The selected components are `icomp` through `icomp + ncomp - 1`.
 *
 * This is typically called for the sources of the Maxwell equations (J/rho)
 * after deposition from the macroparticles.
 *
 * This updates both the *valid* cells and *guard* cells in the selected
 * component range.
 *
 * \param[in] icomp First component of `mf` to sum.
 * \param[in] ncomp Number of components to process.
 */
void
WarpXSumGuardCells(amrex::MultiFab& mf, const amrex::Periodicity& period,
                   const amrex::IntVect& src_ngrow, int icomp, int ncomp);

/** \brief Sum overlapping cell values in `src` independently for each
 * selected component and store the results in `dst`.
 *
 * Source components 0 through `ncomp - 1` map to destination components
 * `icomp` through `icomp + ncomp - 1`, respectively.
 *
 * This is typically called for the sources of the Maxwell equations (J/rho)
 * after deposition from the macroparticles + filtering.
 *
 * This overwrites both the *valid* cells and *guard* cells of `dst` in the
 * selected component range.
 *
 * \param[in] icomp First component of `dst` in which to store the sums.
 * \param[in] ncomp Number of components to process.
 */
void
WarpXSumGuardCells(amrex::MultiFab& dst, amrex::MultiFab& src,
                   const amrex::Periodicity& period,
                   const amrex::IntVect& src_ngrow, int icomp, int ncomp);

#endif // WARPX_SUM_GUARD_CELLS_H_
