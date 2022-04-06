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

    h_quad.reset(new HardEdgedQuadrupole());
    if (h_quad->nelements > 0) {
#ifdef AMREX_USE_GPU
        d_quad = static_cast<HardEdgedQuadrupole*>
            (amrex::The_Arena()->alloc(sizeof(HardEdgedQuadrupole)));
        amrex::Gpu::htod_memcpy_async(d_quad, h_quad.get(), sizeof(HardEdgedQuadrupole));
#else
        d_quad = h_quad.get();
#endif

    }

}

