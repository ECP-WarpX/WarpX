/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "LatticeElementFinder.H"
#include "LatticeElements/HardEdgedQuadrupole.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

using namespace amrex::literals;

LatticeElementFinder::LatticeElementFinder (WarpXParIter const& a_pti)
{

    // The lattice is assumed to extend in the z-direction
    // Get the number of cells where indices will be setup
    amrex::Box box = a_pti.tilebox();
    nz = box.size()[WARPX_ZINDEX] - 1;

    int const lev = 0;
    zmin = WarpX::LowerCorner(box, lev, 0._rt)[2];
    dz = WarpX::CellSize(lev)[2];

    Setup_Quad();

}

void
LatticeElementFinder::Setup_Quad()
{
    // Grab pointers to lattice elements
    auto& warpx = WarpX::GetInstance();
    d_quad = warpx.m_accelerator_lattice->d_quad;

    if (d_quad) {
        // Allocate the space for the indices.
        // Note that this uses nz since the information is saved per cell
        // The "-1" is the flag that no elements at that location
        quad_indices.resize(nz, -1);

        setup_lattice_indices(d_quad->d_zs, d_quad->d_ze, quad_indices);
    }

}
