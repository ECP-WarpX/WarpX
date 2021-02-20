/* Copyright 2021 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARIES_H_
#define PARTICLEBOUNDARIES_H_

#include "Utils/WarpXUtil.H"

enum struct ParticleBC { none=0, absorbing, reflecting };

/**
 *  This allows a different boundary condition on each edge of the grid.
 */
struct ParticleBoundaries
{

    ParticleBoundaries () noexcept
    {
        SetAll("none");
    }

    std::map<std::string, ParticleBC> m_BCname_to_enum = {
        {"none", ParticleBC::none},
        {"absorbing", ParticleBC::absorbing},
        {"reflecting", ParticleBC::reflecting}
    };

    void
    SetAll (std::string bc_name) {
        auto it = m_BCname_to_enum.find(bc_name);
        WarpXUtilMsg::AlwaysAssert(it != m_BCname_to_enum.end(),
                                   "ERROR: unknown particle boundary type '" + bc_name + "'");
        ParticleBC bc = it->second;
        xmin_bc = bc;
        xmax_bc = bc;
        ymin_bc = bc;
        ymax_bc = bc;
        zmin_bc = bc;
        zmax_bc = bc;
    }

    void
    SetBoundsX(std::vector<std::string> bc_names) { SetBounds(bc_names, xmin_bc, xmax_bc); }
    void
    SetBoundsY(std::vector<std::string> bc_names) { SetBounds(bc_names, ymin_bc, ymax_bc); }
    void
    SetBoundsZ(std::vector<std::string> bc_names) { SetBounds(bc_names, zmin_bc, zmax_bc); }

    bool
    AllNone () {
        return (xmin_bc == ParticleBC::none && xmax_bc == ParticleBC::none
             && ymin_bc == ParticleBC::none && ymax_bc == ParticleBC::none
             && zmin_bc == ParticleBC::none && zmax_bc == ParticleBC::none);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void
    apply_boundary_x (amrex::ParticleReal& x, amrex::ParticleReal& v, bool& particle_lost,
                      amrex::Real xmin, amrex::Real xmax, ParticleBoundaries const& boundaries)
    {
        apply_boundary(x, v, particle_lost, xmin, xmax, boundaries.xmin_bc, boundaries.xmax_bc);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void
    apply_boundary_y (amrex::ParticleReal& y, amrex::ParticleReal& v, bool& particle_lost,
                      amrex::Real ymin, amrex::Real ymax, ParticleBoundaries const& boundaries)
    {
        apply_boundary(y, v, particle_lost, ymin, ymax, boundaries.ymin_bc, boundaries.ymax_bc);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void
    apply_boundary_z (amrex::ParticleReal& z, amrex::ParticleReal& v, bool& particle_lost,
                      amrex::Real zmin, amrex::Real zmax, ParticleBoundaries const& boundaries)
    {
        apply_boundary(z, v, particle_lost, zmin, zmax, boundaries.zmin_bc, boundaries.zmax_bc);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void
    apply_boundary (amrex::ParticleReal& x, amrex::ParticleReal& v, bool& particle_lost,
                    amrex::Real xmin, amrex::Real xmax,
                    ParticleBC xmin_bc, ParticleBC xmax_bc)
    {
        if (x < xmin) {
            if (xmin_bc == ParticleBC::absorbing) {
                particle_lost = true;
            }
            else if (xmin_bc == ParticleBC::reflecting) {
                x = 2*xmin - x;
                v = -v;
            }
        }
        else if (x > xmax) {
            if (xmax_bc == ParticleBC::absorbing) {
                particle_lost = true;
            }
            else if (xmax_bc == ParticleBC::reflecting) {
                x = 2*xmax - x;
                v = -v;
            }
        }
    }

private:

    void
    SetBounds(std::vector<std::string> bc_names, ParticleBC& min_bc, ParticleBC& max_bc) {
        if (bc_names.size() == 0) return;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bc_names.size() <= 2,
                                         "ERROR: Too many particle bounderies specified");

        auto it = m_BCname_to_enum.find(bc_names[0]);
        WarpXUtilMsg::AlwaysAssert(it != m_BCname_to_enum.end(),
                                   "ERROR: unknown particle boundary type '" + bc_names[0] + "'");
        min_bc = it->second;
        if (bc_names.size() == 1) {
            max_bc = min_bc;
        } else {
            auto it2 = m_BCname_to_enum.find(bc_names[1]);
            WarpXUtilMsg::AlwaysAssert(it2 != m_BCname_to_enum.end(),
                                       "ERROR: unknown particle boundary type '" + bc_names[1] + "'");
            max_bc = it2->second;
        }
    }

    ParticleBC xmin_bc;
    ParticleBC xmax_bc;
    ParticleBC ymin_bc;
    ParticleBC ymax_bc;
    ParticleBC zmin_bc;
    ParticleBC zmax_bc;

};
#endif /*PARTICLEBOUNDARIES_H_*/
