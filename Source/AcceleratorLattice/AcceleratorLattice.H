/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_
#define WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_

#include "LatticeElements/LatticeElement.H"
#include "LatticeElements/HardEdgedQuadrupole.H"

#include <AMReX_Vector.H>

/* \brief AcceleratorLattice contains the lattice elements
 */
class AcceleratorLattice
{
public:
    AcceleratorLattice ();

    amrex::Vector< LatticeElement* > all_elements;

private:

    HardEdgedQuadrupole quad;

};

#endif // WARPX_ACCELERATORLATTICE_ACCELERATORLATTICE_H_
