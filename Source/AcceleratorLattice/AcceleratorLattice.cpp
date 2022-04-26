/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "AcceleratorLattice.H"
#include "LatticeElements/HardEdgedQuadrupole.H"

AcceleratorLattice::AcceleratorLattice ()
{

    h_quad = std::make_unique<HardEdgedQuadrupole>();
    if (h_quad->nelements > 0) {
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
        UpdateElementFinder();
    }
}

void
AcceleratorLattice::UpdateElementFinder ()
{
    if (m_lattice_defined) {
        for (amrex::MFIter mfi(*m_element_finder); mfi.isValid(); ++mfi)
        {
            (*m_element_finder)[mfi].UpdateIndices(*this);
        }
    }
}

LatticeElementFinderDevice
AcceleratorLattice::Get_Finder_Device(WarpXParIter const& a_pti, int const a_offset) const
{
    LatticeElementFinder & finder = (*m_element_finder)[a_pti];
    return finder.Get_Finder_Device(a_pti, a_offset, *this);
}
