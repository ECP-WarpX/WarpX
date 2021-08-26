/* Copyright 2021 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARIES_K_H_
#define PARTICLEBOUNDARIES_K_H_

#include "ParticleBoundaries.H"

#include <AMReX_AmrCore.H>

namespace ApplyParticleBoundaries {

    /* \brief Applies the boundary condition on a specific axis
     *        This is called by apply_boundaries.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    apply_boundary (amrex::ParticleReal& x, amrex::Real xmin, amrex::Real xmax,
                    bool& change_sign_ux, bool& particle_lost,
                    ParticleBoundaryType xmin_bc, ParticleBoundaryType xmax_bc)
    {
        if (x < xmin) {
            if (xmin_bc == ParticleBoundaryType::Absorbing || xmin_bc == ParticleBoundaryType::Open) {
                particle_lost = true;
            }
            else if (xmin_bc == ParticleBoundaryType::Reflecting) {
                x = 2*xmin - x;
                change_sign_ux = true;
            }
        }
        else if (x > xmax) {
            if (xmax_bc == ParticleBoundaryType::Absorbing || xmin_bc == ParticleBoundaryType::Open) {
                particle_lost = true;
            }
            else if (xmax_bc == ParticleBoundaryType::Reflecting) {
                x = 2*xmax - x;
                change_sign_ux = true;
            }
        }
    }

    /* \brief Applies absorbing or reflecting boundary condition to the input particles, along all axis.
     *        For reflecting boundaries, the position of the particle is changed appropriately and
     *        the sign of the velocity is changed (depending on the reflect_all_velocities flag).
     *        For absorbing, a flag is set whether the particle has been lost (it is up to the calling
     *        code to take appropriate action to remove any lost particles).
     *        Note that periodic boundaries are handled in AMReX code.
     *
     * \param x, xmin, xmax: particle x position, location of x boundary
     * \param y, ymin, ymax: particle y position, location of y boundary (3D only)
     * \param z, zmin, zmax: particle z position, location of z boundary
     * \param ux, uy, uz: particle momenta
     * \param particle_lost: output, flags whether the particle was lost
     * \param boundaries: object with boundary condition settings
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    apply_boundaries (amrex::ParticleReal& x, amrex::Real xmin, amrex::Real xmax,
#ifdef WARPX_DIM_3D
                      amrex::ParticleReal& y, amrex::Real ymin, amrex::Real ymax,
#endif
                      amrex::ParticleReal& z, amrex::Real zmin, amrex::Real zmax,
                      amrex::ParticleReal& ux, amrex::ParticleReal& uy, amrex::ParticleReal& uz,
                      bool& particle_lost,
                      ParticleBoundaries const& boundaries)
    {
        bool change_sign_ux = false;
        bool change_sign_uy = false;
        bool change_sign_uz = false;

        apply_boundary(x, xmin, xmax, change_sign_ux, particle_lost, boundaries.xmin_bc, boundaries.xmax_bc);
#ifdef WARPX_DIM_3D
        apply_boundary(y, ymin, ymax, change_sign_uy, particle_lost, boundaries.ymin_bc, boundaries.ymax_bc);
#endif
        apply_boundary(z, zmin, zmax, change_sign_uz, particle_lost, boundaries.zmin_bc, boundaries.zmax_bc);

        if (boundaries.reflect_all_velocities && (change_sign_ux | change_sign_uy | change_sign_uz)) {
            change_sign_ux = true;
            change_sign_uy = true;
            change_sign_uz = true;
        }
        if (change_sign_ux) ux = -ux;
        if (change_sign_uy) uy = -uy;
        if (change_sign_uz) uz = -uz;
    }

}
#endif
