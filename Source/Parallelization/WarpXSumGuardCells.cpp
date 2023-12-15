/* Copyright 2019 Maxence Thevenet, Remi Lehe, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpXSumGuardCells.H"

#include "Utils/WarpXAlgorithmSelection.H"

#include "WarpX.H"

#include <ablastr/utils/Communication.H>

void
WarpXSumGuardCells(amrex::MultiFab& mf, const amrex::Periodicity& period,
                   const amrex::IntVect& src_ngrow,
                   const int icomp, const int ncomp)
{
    amrex::IntVect n_updated_guards;

    // Update both valid cells and guard cells
    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD) {
        n_updated_guards = mf.nGrowVect();
    } else {  // Update only the valid cells
        n_updated_guards = amrex::IntVect::TheZeroVector();
    }
    ablastr::utils::communication::SumBoundary(mf, icomp, ncomp, src_ngrow, n_updated_guards, WarpX::do_single_precision_comms, period);
}


void
WarpXSumGuardCells(amrex::MultiFab& dst, amrex::MultiFab& src,
                   const amrex::Periodicity& period,
                   const amrex::IntVect& src_ngrow,
                   const int icomp, const int ncomp)
{
    amrex::IntVect n_updated_guards;

    // Update both valid cells and guard cells
    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD) {
        n_updated_guards = dst.nGrowVect();
    } else {  // Update only the valid cells
        n_updated_guards = amrex::IntVect::TheZeroVector();
    }

    dst.setVal(0., icomp, ncomp, n_updated_guards);
//    ablastr::utils::communication::ParallelAdd(dst, src, 0, icomp, ncomp, src_ngrow, n_updated_guards, period);
    dst.ParallelAdd(src, 0, icomp, ncomp, src_ngrow, n_updated_guards, period);
}
