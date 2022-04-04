/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "LatticeElementFinder.H"
#include "LatticeElements/LatticeElement.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

using namespace amrex::literals;

LatticeElementFinder::LatticeElementFinder (WarpXParIter const& a_pti)
{

    auto& warpx = WarpX::GetInstance();

    // Grab a shared pointer
    accelerator_lattice = warpx.m_accelerator_lattice;

    nlattices = static_cast<int>(accelerator_lattice->all_elements.size());

    // If no lattice element types have been setup, then do nothing
    if (nlattices == 0) return;

    // The lattice is assumed to extend in the z-direction
    // Get the number of cells where indices will be setup
    amrex::Box box = a_pti.tilebox();
    nz = box.size()[WARPX_ZINDEX] - 1;

    // Allocate the space for the indices.
    // There will be a set of nz indices for each element type.
    lattice_indices.resize(nlattices);
    for (int itype=0 ; itype < nlattices ; itype++) {
        // Note that this uses nz since the information is saved per cell
        // The "-1" is the flag that no elements at that location
        lattice_indices[itype].resize(nz, -1);
    }

    int const lev = 0;
    zmin = WarpX::LowerCorner(box, lev, 0._rt)[2];
    dz = WarpX::CellSize(lev)[2];

    setup_lattice_indices();

}
