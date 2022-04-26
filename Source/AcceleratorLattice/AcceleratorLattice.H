/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_
#define WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_

#include "LatticeElementFinder.H"
#include "LatticeElements/HardEdgedQuadrupole.H"

/* \brief AcceleratorLattice contains the lattice elements
 */
class AcceleratorLattice
{
public:

    AcceleratorLattice ();

    int m_level;
    bool m_lattice_defined = false;

    void InitElementFinder (int const lev, amrex::BoxArray const & ba, amrex::DistributionMapping const & dm);
    void UpdateElementFinder (int const lev);

    std::unique_ptr<amrex::LayoutData<LatticeElementFinder>> m_element_finder;
    LatticeElementFinderDevice GetFinderDeviceInstance (WarpXParIter const& a_pti, int const a_offset) const;

    std::unique_ptr<HardEdgedQuadrupole> h_quad;

};

#endif // WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_
