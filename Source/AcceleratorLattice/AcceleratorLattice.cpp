/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "AcceleratorLattice.H"

AcceleratorLattice::AcceleratorLattice ()
{

    if (quad.isDefined()) {
        all_elements.push_back(&quad);
    }

}

