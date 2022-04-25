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

using namespace amrex::literals;

void
LatticeElementFinder::InitElementFinder (int const lev, amrex::MFIter const& a_mfi,
                                         AcceleratorLattice const& accelerator_lattice)
{

    // The lattice is assumed to extend in the z-direction
    // Get the number of nodes where indices will be setup
    amrex::Box box = a_mfi.tilebox();
    m_nz = box.size()[WARPX_ZINDEX];

    m_zmin = WarpX::LowerCorner(box, lev, 0._rt)[2];
    m_dz = WarpX::CellSize(lev)[2];

    AllocateQuadIndex(accelerator_lattice);

    UpdateIndices(accelerator_lattice);

}

void
LatticeElementFinder::UpdateIndices (AcceleratorLattice const& accelerator_lattice)
{
    UpdateQuadIndex(accelerator_lattice);
}

void
LatticeElementFinder::AllocateQuadIndex(AcceleratorLattice const& accelerator_lattice)
{
    // Grab pointers to lattice elements
    HardEdgedQuadrupole const *h_quad = accelerator_lattice.h_quad.get();

    if (h_quad) {
        // Allocate the space for the indices.
        // Note that this uses m_nz since the information is saved per cell
        // The "-1" is the flag that no elements at that location
        d_quad_indices.resize(m_nz, -1);
    }
}

void
LatticeElementFinder::UpdateQuadIndex(AcceleratorLattice const& accelerator_lattice)
{
    // Grab pointers to lattice elements
    HardEdgedQuadrupole const *h_quad = accelerator_lattice.h_quad.get();

    if (h_quad) {
        setup_lattice_indices(h_quad->d_zs, h_quad->d_ze, d_quad_indices);
    }
}

LatticeElementFinderDevice
LatticeElementFinder::Get_Finder_Device (WarpXParIter const& a_pti, int const a_offset,
                                         AcceleratorLattice const& accelerator_lattice)
{
    LatticeElementFinderDevice result;
    result.InitLatticeElementFinderDevice(a_pti, a_offset, accelerator_lattice, *this);
    return result;
}


void
LatticeElementFinderDevice::InitLatticeElementFinderDevice (WarpXParIter const& a_pti, int const a_offset,
                                                            AcceleratorLattice const& accelerator_lattice,
                                                            LatticeElementFinder const & h_finder)
{

    auto& warpx = WarpX::GetInstance();

    int const lev = a_pti.GetLevel();

    m_get_position = GetParticlePosition(a_pti, a_offset);
    auto& attribs = a_pti.GetAttribs();
    m_ux = attribs[PIdx::ux].dataPtr() + a_offset;
    m_uy = attribs[PIdx::uy].dataPtr() + a_offset;
    m_uz = attribs[PIdx::uz].dataPtr() + a_offset;
    m_dt = warpx.getdt(lev);

    m_zmin = h_finder.m_zmin;
    m_dz = h_finder.m_dz;

    HardEdgedQuadrupole const *h_quad = accelerator_lattice.h_quad.get();
    if (h_quad) {
        d_quad = h_quad->GetDevice();
        d_quad_indices_arr = h_finder.d_quad_indices.data();
    }

}

