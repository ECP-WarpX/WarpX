/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_

#include "LatticeElements/HardEdgedQuadrupole.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>

class AcceleratorLattice;
struct LatticeElementFinderDevice;

// The top layer needs to be saved in the AcceleratorLattice class in a Layout.
// This gives a place to store the index array device pointers.

struct LatticeElementFinder
{

    void InitElementFinder (int const lev, amrex::MFIter const& a_mfi,
                            AcceleratorLattice const& accelerator_lattice);
    void AllocateIndices (AcceleratorLattice const& accelerator_lattice);
    void UpdateIndices (int const lev, amrex::MFIter const& a_mfi,
                        AcceleratorLattice const& accelerator_lattice);

    int m_nz;
    amrex::Real m_zmin;
    amrex::Real m_dz;

    LatticeElementFinderDevice GetFinderDeviceInstance (WarpXParIter const& a_pti, int const a_offset,
                                                        AcceleratorLattice const& accelerator_lattice);

    amrex::Gpu::DeviceVector<int> d_quad_indices;

    void setup_lattice_indices (amrex::Gpu::DeviceVector<amrex::Real> const & zs,
                                amrex::Gpu::DeviceVector<amrex::Real> const & ze,
                                amrex::Gpu::DeviceVector<int> & indices)
    {

        using namespace amrex::literals;

        int nelements = static_cast<int>(zs.size());
        amrex::Real const * zs_arr = zs.data();
        amrex::Real const * ze_arr = ze.data();
        int * indices_arr = indices.data();

        amrex::Real zmin = m_zmin;
        amrex::Real dz = m_dz;

        amrex::ParallelFor( m_nz,
            [=] AMREX_GPU_DEVICE (int iz) {

                // Get the location of the grid node
                amrex::Real const z_node = zmin + iz*dz;

                // Check if any elements overlap the grid cell, and if so, store its index
                // For now, this assumes that there is no overlap among elements of the same type
                for (int ie = 0 ; ie < nelements ; ie++) {
                    amrex::Real zcenter_left, zcenter_right;
                    if (ie == 0) {
                        zcenter_left = std::numeric_limits<amrex::Real>::lowest();
                    } else {
                        zcenter_left = 0.5_rt*(ze_arr[ie-1] + zs_arr[ie]);
                    }
                    if (ie < nelements - 1) {
                        zcenter_right = 0.5_rt*(ze_arr[ie] + zs_arr[ie+1]);
                    } else {
                        zcenter_right = std::numeric_limits<amrex::Real>::max();
                    }

                    if (zcenter_left <= z_node && z_node < zcenter_right) {
                        indices_arr[iz] = ie;
                    }

                }
            }
        );
    }

};

struct LatticeElementFinderDevice
{

    void
    InitLatticeElementFinderDevice (WarpXParIter const& a_pti, int const a_offset,
                                    AcceleratorLattice const& accelerator_lattice,
                                    LatticeElementFinder const & h_finder);

    amrex::Real m_gamma_boost;
    amrex::Real m_uz_boost;
    amrex::Real m_time;

    amrex::Real m_zmin;
    amrex::Real m_dz;
    amrex::Real m_dt;

    GetParticlePosition m_get_position;
    const amrex::ParticleReal* AMREX_RESTRICT m_ux = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uy = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uz = nullptr;

    HardEdgedQuadrupoleDevice d_quad;
    int const* d_quad_indices_arr = nullptr;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (const long i,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& field_Ez,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& field_Bz) const noexcept
    {

        using namespace amrex::literals;

        amrex::ParticleReal x, y, z;
        m_get_position(i, x, y, z);

        constexpr amrex::Real inv_c2 = 1._rt/(PhysConst::c*PhysConst::c);

        if (m_gamma_boost > 1._rt) {
            z = m_gamma_boost*z + m_uz_boost*m_time;
        }

        const int iz = static_cast<int>((z - m_zmin)/m_dz);

        amrex::Real Ex_sum = 0._rt;
        amrex::Real Ey_sum = 0._rt;
        amrex::Real Ez_sum = 0._rt;
        amrex::Real Bx_sum = 0._rt;
        amrex::Real By_sum = 0._rt;
        amrex::Real Bz_sum = 0._rt;

        if (d_quad.nelements > 0) {
            if (d_quad_indices_arr[iz] > -1) {
                int ielement = d_quad_indices_arr[iz];
                amrex::Real Ex, Ey, Bx, By;
                d_quad.get_field(ielement, x, y, z, m_ux[i], m_uy[i], m_uz[i], m_dt, Ex, Ey, Bx, By);
                Ex_sum += Ex;
                Ey_sum += Ey;
                Bx_sum += Bx;
                By_sum += By;
            }
        }

        if (m_gamma_boost > 1._rt) {
            // Transform the fields to the boosted frame
            const amrex::Real Ex_boost = m_gamma_boost*Ex_sum - m_uz_boost*By_sum;
            const amrex::Real Ey_boost = m_gamma_boost*Ey_sum + m_uz_boost*Bx_sum;
            const amrex::Real Bx_boost = m_gamma_boost*Bx_sum + m_uz_boost*Ey_sum*inv_c2;
            const amrex::Real By_boost = m_gamma_boost*By_sum - m_uz_boost*Ex_sum*inv_c2;
            Ex_sum = Ex_boost;
            Ey_sum = Ey_boost;
            Bx_sum = Bx_boost;
            By_sum = By_boost;
        }

        field_Ex += Ex_sum;
        field_Ey += Ey_sum;
        field_Ez += Ez_sum;
        field_Bx += Bx_sum;
        field_By += By_sum;
        field_Bz += Bz_sum;

    }

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENTFINDER_H_
