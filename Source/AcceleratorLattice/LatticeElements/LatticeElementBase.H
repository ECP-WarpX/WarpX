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

    /**
     * \brief LatticeElementBase constructor
     * This gets the input z extent of the element and check for its correctness
     *
     * @param[in] element_name the name of the element type
     */
    LatticeElementBase (std::string const& element_name);

    ~LatticeElementBase () = default;

    /**
     * \brief Check the correctness of the user input z extents of the elements
     *
     * @param[in] zs the z starts of the elements
     * @param[in] ze the z ends of the elements
     */
    void
    CheckElementCorrectness (amrex::Vector<amrex::Real> const & zs,
                             amrex::Vector<amrex::Real> const & ze);

    std::string m_element_name;

    int nelements = 0;

    // The device level copy of the data.
    // Note that the data is not needed at the host level.
    amrex::Gpu::DeviceVector<amrex::Real> d_zs;
    amrex::Gpu::DeviceVector<amrex::Real> d_ze;

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTBASE_H_
