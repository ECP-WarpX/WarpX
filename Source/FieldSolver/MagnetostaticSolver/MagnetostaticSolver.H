/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: S. Eric Clark (LLNL)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MAGNETOSTATICSOLVER_H_
#define WARPX_MAGNETOSTATICSOLVER_H_

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_REAL.H>


namespace MagnetostaticSolver {

    /** Boundary Handler for the Vector Potential Poisson Solver
     * This only will handle homogeneous Dirichlet boundary conditions on
     * embedded boundaries, and homogeneous dirichlet/Neumann or periodic boundary conditions
     */
    class VectorPoissonBoundaryHandler {
    public:
        amrex::Array<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>, 3> lobc, hibc;
        bool bcs_set = false;
        std::array<std::array<bool, AMREX_SPACEDIM * 2>, 3> dirichlet_flag;
        bool has_non_periodic = false;

        void defineVectorPotentialBCs ();
    };
} // namespace MagnetostaticSolver

#endif //WARPX_MAGNETOSTATICSOLVER_H_
