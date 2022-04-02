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

class LatticeElementFinder
{
public:

    LatticeElementFinder () = default;

    LatticeElementFinder (WarpXParIter const& a_pti, int const a_offset);

    virtual ~LatticeElementFinder () = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (long i,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& field_Ez,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& field_Bz) const noexcept
    {

        amrex::ParticleReal x, y, z;
        m_get_position(i, x, y, z);

        const int iz = static_cast<int>((z - zmin)/dz);

        // Loop over the element type and get the field from any elements that
        // the particle is in.
        int nlattices = static_cast<int>(accelerator_lattice->all_elements.size());
        for (int itype = 0 ; itype < nlattices ; itype++) {
            if (lattice_indices[itype][iz] > -1) {
                amrex::ParticleReal Ex, Ey, Ez, Bx, By, Bz;
                int const ielement = lattice_indices[itype][iz];
                auto element = accelerator_lattice->all_elements[itype];
                element->get_field(ielement, x, y, z, Ex, Ey, Ez, Bx, By, Bz);
                field_Ex += Ex;
                field_Ey += Ey;
                field_Ez += Ez;
                field_Bx += Bx;
                field_By += By;
                field_Bz += Bz;
            }
        }
    }

private:

    amrex::Real zmin;
    amrex::Real dz;

    // Outer vector is types of elements
    // Inner vector is z cells
    amrex::Vector< amrex::Vector <int> > lattice_indices;

    std::shared_ptr<AcceleratorLattice> accelerator_lattice;

    GetParticlePosition m_get_position;



};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
