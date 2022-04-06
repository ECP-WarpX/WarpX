/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_

#include "AcceleratorLattice.H"
#include "LatticeElements/HardEdgedQuadrupole.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>

#include <string>

struct LatticeElementFinder
{

    LatticeElementFinder () = default;

    LatticeElementFinder (WarpXParIter const& a_pti);

    virtual ~LatticeElementFinder () = default;

    int nz;
    amrex::Real zmin;
    amrex::Real dz;

    HardEdgedQuadrupole* d_quad;
    amrex::Gpu::DeviceVector<int> quad_indices;
    void Setup_Quad();

    void setup_lattice_indices(amrex::Gpu::DeviceVector<amrex::Real> & zs,
                               amrex::Gpu::DeviceVector<amrex::Real> & ze,
                               amrex::Gpu::DeviceVector<int> & indices)
    {

        amrex::Real const * zs_arr = zs.data();
        amrex::Real const * ze_arr = ze.data();

        int * indices_arr = indices.data();
        int nelements = static_cast<int>(zs.size());

        amrex::ParallelFor( nz,
            [=] AMREX_GPU_DEVICE (int iz) {

                // Get the location of the grid cell
                amrex::Real const z_lower = zmin + iz*dz;
                amrex::Real const z_upper = zmin + iz*dz + dz;

                // Check if any elements overlap the grid cell, and if so, store its index
                // For now, this assumes that there is no overlap among elements of the same type
                for (int ie = 0 ; ie < nelements ; ie++) {
                    if (zs_arr[ie] <= z_upper && z_lower <= ze_arr[ie]) {
                        indices_arr[iz] = ie;
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
                      amrex::ParticleReal& /*field_Ez*/,
                      amrex::ParticleReal& /*field_Bx*/,
                      amrex::ParticleReal& /*field_By*/,
                      amrex::ParticleReal& /*field_Bz*/) const noexcept
    {

        amrex::ParticleReal x, y, z;
        get_position(i, x, y, z);

        const int iz = static_cast<int>((z - zmin)/dz);

        if (d_quad) {
            if (quad_indices[iz] > -1) {
                int ielement = quad_indices[iz];
                amrex::Real Ex, Ey;
                d_quad->get_field(ielement, x, y, z, ux[i], uy[i], uz[i], dt, Ex, Ey);
                field_Ex += Ex;
                field_Ey += Ey;
            }
        }

    }

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
