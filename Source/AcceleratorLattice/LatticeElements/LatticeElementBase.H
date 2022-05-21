/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTBASE_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTBASE_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>

#include <string>

struct LatticeElementBase
{

    LatticeElementBase (std::string const& element_name);

    ~LatticeElementBase () = default;

    void
    CheckElementCorrectness (amrex::Vector<amrex::Real> const & zs,
                             amrex::Vector<amrex::Real> const & ze);

    std::string m_element_name;

    int nelements = 0;

    amrex::Gpu::DeviceVector<amrex::Real> d_zs;
    amrex::Gpu::DeviceVector<amrex::Real> d_ze;

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTBASE_H_
