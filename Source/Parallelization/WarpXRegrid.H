/* Copyright 2022 Andrew Myers, Axel Huebl, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARALLELIZATION_REGRID_H
#define WARPX_PARALLELIZATION_REGRID_H

#include <AMReX.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <memory>
#include <utility>
#include <vector>


/** ...
 *
 * @tparam T_MultiFabType ..
 * @param mf ..
 * @param dm ..
 * @param redistribute ..
 */
template< typename T_MultiFabType >
void
RemakeMultiFab (
    T_MultiFabType & mf,
    amrex::DistributionMapping const& dm,
    const bool redistribute)
{
    const amrex::IntVect& ng = mf.nGrowVect();
    auto pmf = T_MultiFabType(mf.boxArray(), dm, mf.nComp(), ng);
    if (redistribute) pmf.Redistribute(mf, 0, 0, mf.nComp(), ng);
    mf = std::move(pmf);
}

#endif // WARPX_PARALLELIZATION_REGRID_H
