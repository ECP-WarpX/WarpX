/* Copyright 2022 S. Eric Clark, LLNL
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef MAGNETOSTATICSOLVER_H_
#define MAGNETOSTATICSOLVER_H_

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

    /** use amrex to directly calculate the magnetic field since with EB's the
     *
     * simple finite difference scheme in WarpX::computeE sometimes fails
     */
    class EBCalcBfromVectorPotentialPerLevel {
      private:
        const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>>& m_b_field;
        const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>>& m_grad_buf_e_stag;
        const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>>& m_grad_buf_b_stag;

      public:
        EBCalcBfromVectorPotentialPerLevel(const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>>& b_field,
                                           const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>>& grad_buf_e_stag,
                                           const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>>& grad_buf_b_stag)
                : m_b_field(b_field),
                  m_grad_buf_e_stag(grad_buf_e_stag),
                  m_grad_buf_b_stag(grad_buf_b_stag)
        {}

        void operator()(amrex::Array<std::unique_ptr<amrex::MLMG>,3> & mlmg, int const lev);

        // Function to perform interpolation from cell edges to cell faces
        void doInterp(const std::unique_ptr<amrex::MultiFab> &src, const std::unique_ptr<amrex::MultiFab> &dst);
    };
} // namespace MagnetostaticSolver

#endif //MAGNETOSTATICSOLVER_H_
