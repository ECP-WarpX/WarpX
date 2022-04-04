/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_

#include "AcceleratorLattice.H"
#include "LatticeElements/LatticeElement.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

struct LatticeElementFinder
{

    LatticeElementFinder () = default;

    LatticeElementFinder (WarpXParIter const& a_pti);

    virtual ~LatticeElementFinder () = default;

    AMREX_GPU_HOST_DEVICE
    void setup_lattice_indices()
    {

        amrex::ParallelFor( nz,
            [=] AMREX_GPU_DEVICE (int iz) {

                // Get the location of the grid cell
                amrex::Real const z_lower = zmin + iz*dz;
                amrex::Real const z_upper = zmin + iz*dz + dz;

                // For each grid node along z, find the element of each type that overlaps that grid cell
                // Loop over the element types that have been defined
                for (int itype=0 ; itype < nlattices ; itype++) {
                    auto const element = accelerator_lattice->all_elements[itype];
                    auto const zs = element->get_zstarts().dataPtr();
                    auto const ze = element->get_zends().dataPtr();

                    // Check if any elements overlap the grid cell, and if so, store its index
                    // For now, this assumes that there is no overlap among elements of the same type
                    for (int ie = 0 ; ie < element->nelements ; ie++) {
                        if (zs[ie] <= z_upper && z_lower <= ze[ie]) {
                            lattice_indices[itype][iz] = ie;
                        }
                    }
                }
            }
        );
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (const long i,
                      const GetParticlePosition& get_position,
                      const amrex::ParticleReal* AMREX_RESTRICT ux,
                      const amrex::ParticleReal* AMREX_RESTRICT uy,
                      const amrex::ParticleReal* AMREX_RESTRICT uz,
                      const amrex::Real dt,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& field_Ez,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& field_Bz) const noexcept
    {

        if (nlattices == 0) return;

        amrex::ParticleReal x, y, z;
        get_position(i, x, y, z);

        const int iz = static_cast<int>((z - zmin)/dz);

        // Loop over the element type and get the field from any elements that
        // the particle is in.
        for (int itype = 0 ; itype < nlattices ; itype++) {
            if (lattice_indices[itype][iz] > -1) {
                amrex::ParticleReal Ex, Ey, Ez, Bx, By, Bz;
                int const ielement = lattice_indices[itype][iz];
                auto element = accelerator_lattice->all_elements[itype];
                element->get_field(ielement, x, y, z, ux[i], uy[i], uz[i], dt, Ex, Ey, Ez, Bx, By, Bz);
                field_Ex += Ex;
                field_Ey += Ey;
                field_Ez += Ez;
                field_Bx += Bx;
                field_By += By;
                field_Bz += Bz;
            }
        }
    }

    int nlattices;
    int nz;
    amrex::Real zmin;
    amrex::Real dz;

    // Outer vector is types of elements
    // Inner vector is z cells
    amrex::Vector< amrex::Vector <int> > lattice_indices;

    std::shared_ptr<AcceleratorLattice> accelerator_lattice;

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
