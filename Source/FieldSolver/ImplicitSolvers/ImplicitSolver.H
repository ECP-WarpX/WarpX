/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef Implicit_Solver_H_
#define Implicit_Solver_H_

#include "FieldSolver/ImplicitSolvers/WarpXSolverVec.H"

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

/**
 * \brief Base class for implicit time solvers. The base functions are those
 *  needed by an implicit solver to be used through WarpX and those needed
 *  to use the nonlinear solvers, such as Picard or Newton (i.e., JFNK).
 */

class WarpX;
class ImplicitSolver
{
public:

    ImplicitSolver() = default;

    virtual ~ImplicitSolver() = default;

    // Prohibit Move and Copy operations
    ImplicitSolver(const ImplicitSolver&) = delete;
    ImplicitSolver& operator=(const ImplicitSolver&) = delete;
    ImplicitSolver(ImplicitSolver&&) = delete;
    ImplicitSolver& operator=(ImplicitSolver&&) = delete;

    //
    // the following routines are called by WarpX
    //

    virtual void Define ( WarpX*  a_WarpX ) = 0;

    [[nodiscard]] virtual bool IsDefined () const = 0;

    virtual void PrintParameters () const = 0;

    virtual void Initialize () = 0;

    virtual void GetParticleSolverParams (int&  a_max_particle_iter,
                                          amrex::ParticleReal&  a_particle_tol ) = 0;

    virtual void OneStep ( amrex::Real  a_time,
                           amrex::Real  a_dt,
                           int          a_step ) = 0;

    //
    // the following routines are called by the linear and nonlinear solvers
    //

    virtual void PreRHSOp ( const WarpXSolverVec&  a_E,
                            amrex::Real            a_time,
                            amrex::Real            a_dt,
                            int                    a_nl_iter,
                            bool                   a_from_jacobian ) = 0;

    virtual void ComputeRHS ( WarpXSolverVec&  a_Erhs,
                        const WarpXSolverVec&  a_E,
                              amrex::Real      a_time,
                              amrex::Real      a_dt ) = 0;

private:

};

#endif
