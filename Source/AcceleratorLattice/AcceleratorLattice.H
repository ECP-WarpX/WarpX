/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_
#define WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_

#include "LatticeElements/HardEdgedQuadrupole.H"

/* \brief AcceleratorLattice contains the lattice elements
 */
class AcceleratorLattice
{
public:

    AcceleratorLattice ();

    bool m_lattice_defined = false;

    std::unique_ptr<HardEdgedQuadrupole> h_quad;
    HardEdgedQuadrupole *d_quad = nullptr;

};

#endif // WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_
