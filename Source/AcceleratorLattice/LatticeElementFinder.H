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

    LatticeElementFinder (WarpXParIter const& a_pti, int const a_offset);

    virtual ~LatticeElementFinder () = default;

    bool m_lattice_defined = false;

    int m_nz;
    amrex::Real m_zmin;
    amrex::Real m_dz;
    amrex::Real m_dt;

    GetParticlePosition m_get_position;
    const amrex::ParticleReal* AMREX_RESTRICT m_ux = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uy = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uz = nullptr;

    HardEdgedQuadrupole* d_quad = nullptr;
    amrex::Gpu::DeviceVector<int> quad_indices;
    void Setup_Quad();

    void setup_lattice_indices(amrex::Gpu::DeviceVector<amrex::Real> & zcenters,
                               amrex::Gpu::DeviceVector<int> & indices)
    {

        using namespace amrex::literals;

        int nelements = static_cast<int>(zcenters.size());
        amrex::Real const * zcenters_arr = zcenters.data();
        int * indices_arr = indices.data();

        amrex::ParallelFor( m_nz,
            [=] AMREX_GPU_DEVICE (int iz) {

                // Get the location of the grid node
                amrex::Real const z_node = m_zmin + iz*m_dz;

                // Check if any elements overlap the grid cell, and if so, store its index
                // For now, this assumes that there is no overlap among elements of the same type
                for (int ie = 0 ; ie < nelements ; ie++) {
                    if (zcenters_arr[ie] <= z_node && z_node < zcenters_arr[ie+1]) {
                        indices_arr[iz] = ie;
                    }
                }
            }
        );
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (const long i,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& /*field_Ez*/,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& /*field_Bz*/) const noexcept
    {

        if (!m_lattice_defined) return;

        amrex::ParticleReal x, y, z;
        m_get_position(i, x, y, z);

        const int iz = static_cast<int>((z - m_zmin)/m_dz);

        if (d_quad) {
            if (quad_indices[iz] > -1) {
                int ielement = quad_indices[iz];
                amrex::Real Ex, Ey, Bx, By;
                d_quad->get_field(ielement, x, y, z, m_ux[i], m_uy[i], m_uz[i], m_dt, Ex, Ey, Bx, By);
                field_Ex += Ex;
                field_Ey += Ey;
                field_Bx += Bx;
                field_By += By;
            }
        }

    }

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
