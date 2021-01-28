/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX.H>

#include <array>
#include <cmath>


/**
 * This struct contains only static functions to initialize the stencil coefficients
 * and to compute finite-difference derivatives for the Cartesian Yee algorithm.
 */
struct CylindricalYeeAlgorithm {

    static void InitializeStencilCoefficients (
        std::array<amrex::Real,3>& cell_size,
        amrex::Vector<amrex::Real>& stencil_coefs_r,
        amrex::Vector<amrex::Real>& stencil_coefs_z ) {

        using namespace amrex;
        // Store the inverse cell size along each direction in the coefficients
        stencil_coefs_r.resize(1);
        stencil_coefs_r[0] = 1._rt/cell_size[0];  // 1./dr
        stencil_coefs_z.resize(1);
        stencil_coefs_z[0] = 1._rt/cell_size[2];  // 1./dz
    }

    /** Compute the maximum, CFL-stable timestep
     *
     * Compute the maximum timestep, for which the scheme remains stable
     * under the Courant-Friedrichs-Levy limit.
     */
    static amrex::Real ComputeMaxDt ( amrex::Real const * const dx,
                                      int const n_rz_azimuthal_modes ) {
        using namespace amrex::literals;
        // In the rz case, the Courant limit has been evaluated
        // semi-analytically by R. Lehe, and resulted in the following
        // coefficients.
        std::array< amrex::Real, 6 > const multimode_coeffs = {{ 0.2105, 1.0, 3.5234, 8.5104, 15.5059, 24.5037 }};
        amrex::Real multimode_alpha;
        if (n_rz_azimuthal_modes < 7) {
          // Use the table of the coefficients
          multimode_alpha = multimode_coeffs[n_rz_azimuthal_modes-1];
        } else {
          // Use a realistic extrapolation
          multimode_alpha = (n_rz_azimuthal_modes - 1._rt)*(n_rz_azimuthal_modes - 1._rt) - 0.4_rt;
        }
        amrex::Real delta_t = 1._rt / ( std::sqrt(
                                     (1._rt + multimode_alpha) / (dx[0]*dx[0])
                                    + 1._rt / (dx[1]*dx[1])
                              ) * PhysConst::c );
        return delta_t;
    }

    /** Applies the differential operator `1/r * d(rF)/dr`,
     * where `F` is on a *nodal* grid in `r`
     * and the differential operator is evaluated on a *cell-centered* grid.
     * The input parameter `r` is given at the cell-centered position */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDrr_over_r (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const r, amrex::Real const dr,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        ignore_unused(n_coefs_r);

        Real const inv_dr = coefs_r[0];
        return 1._rt/r * inv_dr*( (r+0.5_rt*dr)*F(i+1,j,k,comp) - (r-0.5_rt*dr)*F(i,j,k,comp) );
    }

    /** Applies the differential operator `1/r * d(rF)/dr`,
     * where `F` is on a *cell-centered* grid in `r`
     * and the differential operator is evaluated on a *nodal* grid.
     * The input parameter `r` is given at the cell-centered position */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDrr_over_r (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const r, amrex::Real const dr,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        ignore_unused(n_coefs_r);

        Real const inv_dr = coefs_r[0];
        return 1._rt/r * inv_dr*( (r+0.5_rt*dr)*F(i,j,k,comp) - (r-0.5_rt*dr)*F(i-1,j,k,comp) );
    }

    /**
     * Perform derivative along r on a cell-centered grid, from a nodal field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDr (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        ignore_unused(n_coefs_r);

        Real const inv_dr = coefs_r[0];
        return inv_dr*( F(i+1,j,k,comp) - F(i,j,k,comp) );
    };

    /**
     * Perform derivative along r on a nodal grid, from a cell-centered field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDr (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        ignore_unused(n_coefs_r);

        Real const inv_dr = coefs_r[0];
        return inv_dr*( F(i,j,k,comp) - F(i-1,j,k,comp) );
    }

    /**
     * Perform derivative along z on a cell-centered grid, from a nodal field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k, int const comp ) {

        amrex::ignore_unused(n_coefs_z);

        amrex::Real const inv_dz = coefs_z[0];
        return inv_dz*( F(i,j+1,k,comp) - F(i,j,k,comp) );
    }

    /**
     * Perform derivative along z on a nodal grid, from a cell-centered field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k, int const comp ) {

        amrex::ignore_unused(n_coefs_z);

        amrex::Real const inv_dz = coefs_z[0];
        return inv_dz*( F(i,j,k,comp) - F(i,j-1,k,comp) );
    }

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_
