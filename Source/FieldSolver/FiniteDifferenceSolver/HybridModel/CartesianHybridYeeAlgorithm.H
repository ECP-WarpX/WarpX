/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_HYBRID_YEE_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_HYBRID_YEE_H_

#include "Utils/CoarsenIO.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#include <array>
#include <cmath>

/**
 * This struct contains only static functions to initialize the stencil coefficients
 * and to compute finite-difference derivatives for the Cartesian Yee algorithm
 * as needed by the hybrid algorithm.
 */
struct CartesianHybridYeeAlgorithm {

    /**
     * \brief Returns maximum number of guard cells required by the field-solve
     */
    // static amrex::IntVect GetMaxGuardCell () {
    //     // The yee solver requires one guard cell in each dimension
    //     return (amrex::IntVect(AMREX_D_DECL(1,1,1)));
    // }

    /**
     * \brief Calculate the x-component of the curl of B (J) with the
     * given staggering.
    */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static amrex::Real Jx (
        amrex::Array4<amrex::Real> const& By,
        amrex::Array4<amrex::Real> const& Bz,
        amrex::GpuArray<int,3> const& By_stag,
        amrex::GpuArray<int,3> const& Bz_stag,
        amrex::Real const * const coefs_y,
        amrex::Real const * const coefs_z,
        amrex::GpuArray<int,3> const& target_stag,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;

        // The "coarsening is just 1 i.e. no coarsening"
        GpuArray<int, 3> const& coarsen = {1,1,1};
        GpuArray<int, 3> dest_stag;

#if (defined WARPX_DIM_1D_Z) || (defined WARPX_DIM_XZ)
        ignore_unused(Bz, Bz_stag, coefs_y);
        auto const dyBz = 0._rt; // y derivative is 0 in 1d and 2d
#else
        Real const inv_dy = coefs_y[0];
        Real dyBz;
        if (target_stag[1] == 1) {
            // if the y component of the target staggering is nodal, the y
            // component of Bz should be cell-centered and a "downward"
            // difference used for the derivative
            dest_stag = {target_stag[0],0,target_stag[2]};
            dyBz = inv_dy*(
                CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i, j, k, ncomp)
                - CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i, j-1, k, ncomp)
            );
        }
        else {
            // if the y component of the target staggering is cell-centered,
            // the y component of Bz should be nodal and a "upward" difference
            // used for the derivative
            dest_stag = {target_stag[0],1,target_stag[2]};
            dyBz = inv_dy*(
                CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i, j+1, k, ncomp)
                - CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i, j, k, ncomp)
            );
        }
#endif

        Real const inv_dz = coefs_z[0];
        Real dzBy;
        dest_stag = target_stag;
        if (target_stag[WARPX_ZINDEX] == 1) {
            // if the z component of the target staggering is nodal, the z
            // component of By should be cell-centered and a "downward"
            // difference used for the derivative
            dest_stag[WARPX_ZINDEX] = 0;
            dzBy = inv_dz*(
                CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j, k, ncomp)
#if (defined WARPX_DIM_1D_Z)
                - CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i-1, j, k, ncomp)
#elif (defined WARPX_DIM_XZ)
                - CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j-1, k, ncomp)
#else
                - CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j, k-1, ncomp)
#endif
            );
        }
        else {
            // if the z component of the target staggering is cell-centered,
            // the z component of By should be nodal and a "upward" difference
            // used for the derivative
            dest_stag[WARPX_ZINDEX] = 1;
            dzBy = inv_dz*(
#if (defined WARPX_DIM_1D_Z)
                CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i+1, j, k, ncomp)
#elif (defined WARPX_DIM_XZ)
                CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j+1, k, ncomp)
#else
                CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j, k+1, ncomp)
#endif
                - CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j, k, ncomp)
            );
        }
        // if (std::abs(dyBz - dzBy) < 1e-6) return 0._rt;
        return (dyBz - dzBy) / PhysConst::mu0;
    }

    /**
     * \brief Calculate the y-component of the curl of B (J) with the
     * given staggering.
    */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static amrex::Real Jy (
        amrex::Array4<amrex::Real> const& Bx,
        amrex::Array4<amrex::Real> const& Bz,
        amrex::GpuArray<int,3> const& Bx_stag,
        amrex::GpuArray<int,3> const& Bz_stag,
        amrex::Real const * const coefs_x,
        amrex::Real const * const coefs_z,
        amrex::GpuArray<int,3> const& target_stag,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;

        // The "coarsening is just 1 i.e. no coarsening"
        GpuArray<int, 3> const& coarsen = {1,1,1};
        GpuArray<int, 3> dest_stag;

#if (defined WARPX_DIM_1D_Z)
        ignore_unused(Bz, Bz_stag, coefs_x);
        auto const dxBz = 0._rt; // x derivative is 0 in 1d
#else
        Real const inv_dx = coefs_x[0];
        Real dxBz;
        dest_stag = target_stag;
        if (target_stag[0] == 1) {
            // if the x component of the target staggering is nodal, the x
            // component of Bz should be cell-centered and a "downward"
            // difference used for the derivative
            dest_stag[0] = 0;
            dxBz = inv_dx*(
                CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i, j, k, ncomp)
                - CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i-1, j, k, ncomp)
            );
        }
        else {
            // if the x component of the target staggering is cell-centered,
            // the x component of Bz should be nodal and a "upward" difference
            // used for the derivative
            dest_stag[0] = 1;
            dxBz = inv_dx*(
                CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i+1, j, k, ncomp)
                - CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i, j, k, ncomp)
            );
        }
#endif

        Real const inv_dz = coefs_z[0];
        Real dzBx;
        dest_stag = target_stag;
        if (target_stag[WARPX_ZINDEX] == 1) {
            // if the z component of the target staggering is nodal, the z
            // component of Bx should be cell-centered and a "downward"
            // difference used for the derivative
            dest_stag[WARPX_ZINDEX] = 0;
            dzBx = inv_dz*(
                CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k, ncomp)
#if (defined WARPX_DIM_1D_Z)
                - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i-1, j, k, ncomp)
#elif (defined WARPX_DIM_XZ)
                - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j-1, k, ncomp)
#else
                - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k-1, ncomp)
#endif
            );
        }
        else {
            // if the z component of the target staggering is cell-centered,
            // the z component of Bx should be nodal and a "upward" difference
            // used for the derivative
            dest_stag[WARPX_ZINDEX] = 1;
            dzBx = inv_dz*(
#if (defined WARPX_DIM_1D_Z)
                CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i+1, j, k, ncomp)
#elif (defined WARPX_DIM_XZ)
                CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j+1, k, ncomp)
#else
                CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k+1, ncomp)
#endif
                - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k, ncomp)
            );
        }
        // if (std::abs(dzBx - dxBz) < 1e-6) return 0._rt;
        return (dzBx - dxBz) / PhysConst::mu0;
    }

    /**
     * \brief Calculate the z-component of the curl of B (J) with the
     * given staggering.
    */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static amrex::Real Jz (
        amrex::Array4<amrex::Real> const& Bx,
        amrex::Array4<amrex::Real> const& By,
        amrex::GpuArray<int,3> const& Bx_stag,
        amrex::GpuArray<int,3> const& By_stag,
        amrex::Real const * const coefs_x,
        amrex::Real const * const coefs_y,
        amrex::GpuArray<int,3> const& target_stag,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;

        // The "coarsening is just 1 i.e. no coarsening"
        GpuArray<int, 3> const& coarsen = {1,1,1};
        GpuArray<int, 3> dest_stag;

#if (defined WARPX_DIM_1D_Z)
        ignore_unused(By, By_stag, coefs_x);
        auto const dxBy = 0._rt; // x derivative is 0 in 1d
#else
        Real const inv_dx = coefs_x[0];
        Real dxBy;
        dest_stag = target_stag;
        if (target_stag[0] == 1) {
            // if the x component of the target staggering is nodal, the x
            // component of By should be cell-centered and a "downward"
            // difference used for the derivative
            dest_stag[0] = 0;
            dxBy = inv_dx*(
                CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j, k, ncomp)
                - CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i-1, j, k, ncomp)
            );
        }
        else {
            // if the x component of the target staggering is cell-centered,
            // the x component of By should be nodal and a "upward" difference
            // used for the derivative
            dest_stag[0] = 1;
            dxBy = inv_dx*(
                CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i+1, j, k, ncomp)
                - CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j, k, ncomp)
            );
        }
#endif

#if (defined WARPX_DIM_1D_Z) || (defined WARPX_DIM_XZ)
        ignore_unused(Bx, Bx_stag, coefs_y);
        auto const dyBx = 0._rt; // the y-derivative is 0 in 1d and 2d
#else
        Real const inv_dy = coefs_y[0];
        Real dyBx;
        dest_stag = target_stag;
        if (target_stag[1] == 1) {
            // if the y component of the target staggering is nodal, the y
            // component of Bx should be cell-centered and a "downward"
            // difference used for the derivative
            dest_stag[1] = 0;
            dyBx = inv_dy*(
                CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k, ncomp)
                - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j-1, k, ncomp)
            );
        }
        else {
            // if the y component of the target staggering is cell-centered,
            // the y component of Bx should be nodal and a "upward" difference
            // used for the derivative
            dest_stag[1] = 1;
            dyBx = inv_dy*(
                CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j+1, k, ncomp)
                - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k, ncomp)
            );
        }
#endif
        // if (std::abs(dxBy - dyBx) < 1e-6) return 0._rt;
        return (dxBy - dyBx) / PhysConst::mu0;
    }

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_HYBRID_YEE_H_