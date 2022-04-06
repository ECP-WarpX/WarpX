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

LatticeElementFinder::LatticeElementFinder (WarpXParIter const& a_pti, int const a_offset)
{

    auto& warpx = WarpX::GetInstance();
    m_lattice_defined = warpx.m_accelerator_lattice->m_lattice_defined;
    if (!m_lattice_defined) return;

    int const lev = a_pti.GetLevel();

    m_get_position = GetParticlePosition(a_pti, a_offset);
    auto& attribs = a_pti.GetAttribs();
    m_ux = attribs[PIdx::ux].dataPtr() + a_offset;
    m_uy = attribs[PIdx::uy].dataPtr() + a_offset;
    m_uz = attribs[PIdx::uz].dataPtr() + a_offset;
    m_dt = warpx.getdt(lev);

    // The lattice is assumed to extend in the z-direction
    // Get the number of nodes where indices will be setup
    amrex::Box box = a_pti.tilebox();
    m_nz = box.size()[WARPX_ZINDEX];

    m_zmin = WarpX::LowerCorner(box, lev, 0._rt)[2];
    m_dz = WarpX::CellSize(lev)[2];

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
        // Note that this uses m_nz since the information is saved per cell
        // The "-1" is the flag that no elements at that location
        quad_indices.resize(m_nz, -1);

        setup_lattice_indices(d_quad->d_zcenters, quad_indices);
    }

}
