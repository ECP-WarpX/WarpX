/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "LatticeElementFinder.H"
#include "LatticeElements/LatticeElement.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

using namespace amrex::literals;

LatticeElementFinder::LatticeElementFinder (WarpXParIter const& a_pti, int const a_offset)
{

    m_get_position = GetParticlePosition(a_pti, a_offset);

    auto& warpx = WarpX::GetInstance();

    // Grab a shared pointer
    accelerator_lattice = warpx.m_accelerator_lattice;

    int nlattices = static_cast<int>(accelerator_lattice->all_elements.size());

    // If no lattice element types have been setup, then do nothing
    if (nlattices == 0) return;

    // The lattice is assumed to extend in the z-direction
    // Get the number of cells where indices will be setup
    amrex::Box box = a_pti.tilebox();
    int nz = box.size()[WARPX_ZINDEX] - 1;

    // Allocate the space for the indices.
    // There will be a set of nz indices for each element type.
    lattice_indices.resize(nlattices);
    for (int itype=0 ; itype < nlattices ; itype++) {
        // Note that this uses nz since the information is saved per cell
        // The "-1" is the flag that no elements at that location
        lattice_indices[itype].resize(nz, -1);
    }

    int const lev = 0;
    zmin = WarpX::LowerCorner(box, lev, 0._rt)[2];
    dz = WarpX::CellSize(lev)[2];

    // For each grid node along z, find the element of each type that overlaps that grid cell
    // Loop over the element types that have been defined
    for (int itype=0 ; itype < nlattices ; itype++) {
        auto const element = accelerator_lattice->all_elements[itype];
        auto const& zs = element->get_zstarts();
        auto const& ze = element->get_zends();

        amrex::ParallelFor( nz,
            [=] AMREX_GPU_DEVICE (int iz) {

                // Get the location of the grid cell
                amrex::Real const z_lower = zmin + iz*dz;
                amrex::Real const z_upper = zmin + iz*dz + dz;

                // Check if any elements overlap the grid cell, and if so, store its index
                // For now, this assumes that there is no overlap among elements of the same type
                for (int ie = 0 ; ie < zs.size() ; ie++) {
                    if (zs[ie] <= z_upper && z_lower <= ze[ie]) {
                        lattice_indices[itype][iz] = ie;
                    }
                }
            }
        );

    }
}
