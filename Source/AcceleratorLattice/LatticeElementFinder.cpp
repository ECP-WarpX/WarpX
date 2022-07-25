/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "LatticeElementFinder.H"
#include "LatticeElements/HardEdgedQuadrupole.H"
#include "LatticeElements/HardEdgedPlasmaLens.H"

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

    m_dz = WarpX::CellSize(lev)[2];

    m_gamma_boost = WarpX::gamma_boost;
    m_uz_boost = std::sqrt(WarpX::gamma_boost*WarpX::gamma_boost - 1._prt)*PhysConst::c;

    AllocateIndices(accelerator_lattice);

    UpdateIndices(lev, a_mfi, accelerator_lattice);

}

void
LatticeElementFinder::AllocateIndices (AcceleratorLattice const& accelerator_lattice)
{
    // Allocate the space for the indices for each element type.
    // Note that this uses m_nz since the information is saved per node.

    if (accelerator_lattice.h_quad.nelements > 0) {
        d_quad_indices.resize(m_nz);
    }

    if (accelerator_lattice.h_plasmalens.nelements > 0) {
        d_plasmalens_indices.resize(m_nz);
    }
}

void
LatticeElementFinder::UpdateIndices (int const lev, amrex::MFIter const& a_mfi,
                                     AcceleratorLattice const& accelerator_lattice)
{
    auto& warpx = WarpX::GetInstance();

    // Update the location of the index grid.
    // Note that the current box is used since the box may have been updated since
    // the initialization in InitElementFinder.
    amrex::Box box = a_mfi.tilebox();
    m_zmin = WarpX::LowerCorner(box, lev, 0._rt)[2];
    m_time = warpx.gett_new(lev);

    if (accelerator_lattice.h_quad.nelements > 0) {
        setup_lattice_indices(accelerator_lattice.h_quad.d_zs,
                              accelerator_lattice.h_quad.d_ze,
                              d_quad_indices);
    }

    if (accelerator_lattice.h_plasmalens.nelements > 0) {
        setup_lattice_indices(accelerator_lattice.h_plasmalens.d_zs,
                              accelerator_lattice.h_plasmalens.d_ze,
                              d_plasmalens_indices);
    }
}

LatticeElementFinderDevice
LatticeElementFinder::GetFinderDeviceInstance (WarpXParIter const& a_pti, int const a_offset,
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

    m_gamma_boost = WarpX::gamma_boost;
    m_uz_boost = std::sqrt(WarpX::gamma_boost*WarpX::gamma_boost - 1._prt)*PhysConst::c;

    m_zmin = h_finder.m_zmin;
    m_dz = h_finder.m_dz;
    m_time = h_finder.m_time;

    if (accelerator_lattice.h_quad.nelements > 0) {
        d_quad = accelerator_lattice.h_quad.GetDeviceInstance();
        d_quad_indices_arr = h_finder.d_quad_indices.data();
    }

    if (accelerator_lattice.h_plasmalens.nelements > 0) {
        d_plasmalens = accelerator_lattice.h_plasmalens.GetDeviceInstance();
        d_plasmalens_indices_arr = h_finder.d_plasmalens_indices.data();
    }

}
