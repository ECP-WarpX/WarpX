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
     * \brief Calculate the z-component of the curl of B with the staggering
     * of the Ex field.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real curlB_z_Ex (
        amrex::Array4<amrex::Real> const& Bx,
        amrex::Array4<amrex::Real> const& By,
        amrex::GpuArray<int,3> const& Bx_stag,
        amrex::GpuArray<int,3> const& By_stag,
        amrex::Real const * const coefs_x,
        amrex::Real const * const coefs_y,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;

        // The "coarsening is just 1 i.e. no coarsening"
        GpuArray<int, 3> const& coarsen = {1,1,1};
        GpuArray<int, 3> dest_stag;

#if (defined WARPX_DIM_1D_Z)
        ignore_unused(Bx, By, Bx_stag, By_stag, coefs_x, coefs_y,
            i, j, k, ncomp);
        return 0._rt; // both x and y derivatives are 0 in 1d
#endif
        // take finite difference derivative along x of By values interpolated to
        // (1, 1, 1) staggering
        Real const inv_dx = coefs_x[0];
        dest_stag = {1,1,1};
        auto dxBy = inv_dx*(
            CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i+1, j, k, ncomp)
            - CoarsenIO::Interp(By, By_stag, dest_stag, coarsen, i, j, k, ncomp)
        );
#if (defined WARPX_DIM_XZ)
        ignore_unused(Bx, Bx_stag, coefs_y);
        return dxBy; // the y-derivative term is zero in 2d
#endif
        Real const inv_dy = coefs_y[0];
        dest_stag = {0,0,1};
        auto dyBx = inv_dy*(
            CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k, ncomp)
            - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j-1, k, ncomp)
        );
        return dxBy - dyBx;
    }


    /**
     * \brief Calculate the y-component of the curl of B with the staggering
     * of the Ex field.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real curlB_y_Ex (
        amrex::Array4<amrex::Real> const& Bx,
        amrex::Array4<amrex::Real> const& Bz,
        amrex::GpuArray<int,3> const& Bx_stag,
        amrex::GpuArray<int,3> const& Bz_stag,
        amrex::Real const * const coefs_x,
        amrex::Real const * const coefs_z,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;

        // The "coarsening is just 1 i.e. no coarsening"
        GpuArray<int, 3> const& coarsen = {1,1,1};
        GpuArray<int, 3> dest_stag;

#if (defined WARPX_DIM_1D_Z)
        ignore_unused(Bz, Bz_stag, coefs_z);
        Real dxBz = 0._rt; // x derivative is 0 in 1d
#else
        // take finite difference derivative along x of Bz values interpolated to
        // (1, 1, 1) staggering
        Real const inv_dx = coefs_x[0];
        dest_stag = {1,1,1};
        auto dxBz = inv_dx*(
            CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i+1, j, k, ncomp)
            - CoarsenIO::Interp(Bz, Bz_stag, dest_stag, coarsen, i, j, k, ncomp)
        );
#endif

        Real const inv_dz = coefs_z[0];
        dest_stag = {0,1,0};
        auto dyBx = inv_dz*(
            CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k, ncomp)
            - CoarsenIO::Interp(Bx, Bx_stag, dest_stag, coarsen, i, j, k-1, ncomp)
        );
        return dxBz - dyBx;
    }

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_HYBRID_YEE_H_