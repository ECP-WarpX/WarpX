/* Copyright 2021 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARIES_H_
#define PARTICLEBOUNDARIES_H_

enum struct ParticleBC { none=0, absorbing, reflecting };

/**
 *  This allows a different boundary condition on each edge of the grid.
 */
struct ParticleBoundaries
{

    ParticleBoundaries () noexcept
    {
        SetAll(ParticleBC::none);
    }

    ParticleBC xmin_bc;
    ParticleBC xmax_bc;
#ifdef WARPX_DIM_3D
    ParticleBC ymin_bc;
    ParticleBC ymax_bc;
#endif
    ParticleBC zmin_bc;
    ParticleBC zmax_bc;

    void
    SetAll (ParticleBC bc) {
        xmin_bc = bc;
        xmax_bc = bc;
#ifdef WARPX_DIM_3D
        ymin_bc = bc;
        ymax_bc = bc;
#endif
        zmin_bc = bc;
        zmax_bc = bc;
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

};
#endif /*PARTICLEBOUNDARIES_H_*/
