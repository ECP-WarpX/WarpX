/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "AcceleratorLattice.H"
#include "LatticeElements/HardEdgedQuadrupole.H"
#include "LatticeElements/HardEdgedPlasmaLens.H"

AcceleratorLattice::AcceleratorLattice ()
{

    /* Get the inputs for and initialize all of the lattice element types */

    h_quad = std::make_unique<HardEdgedQuadrupole>();
    if (h_quad->nelements > 0) {
        m_lattice_defined = true;
    }

    h_plasmalens = std::make_unique<HardEdgedPlasmaLens>();
    if (h_plasmalens->nelements > 0) {
        m_lattice_defined = true;
    }

}

void
AcceleratorLattice::InitElementFinder (int const lev, amrex::BoxArray const & ba, amrex::DistributionMapping const & dm)
{
    if (m_lattice_defined) {
        m_element_finder = std::make_unique<amrex::LayoutData<LatticeElementFinder>>(ba, dm);
        for (amrex::MFIter mfi(*m_element_finder); mfi.isValid(); ++mfi)
        {
            (*m_element_finder)[mfi].InitElementFinder(lev, mfi, *this);
        }
    }
}

void
AcceleratorLattice::UpdateElementFinder (int const lev)
{
    if (m_lattice_defined) {
        for (amrex::MFIter mfi(*m_element_finder); mfi.isValid(); ++mfi)
        {
            (*m_element_finder)[mfi].UpdateIndices(lev, mfi, *this);
        }
    }
}

LatticeElementFinderDevice
AcceleratorLattice::GetFinderDeviceInstance (WarpXParIter const& a_pti, int const a_offset) const
{
    LatticeElementFinder & finder = (*m_element_finder)[a_pti];
    return finder.GetFinderDeviceInstance(a_pti, a_offset, *this);
}
