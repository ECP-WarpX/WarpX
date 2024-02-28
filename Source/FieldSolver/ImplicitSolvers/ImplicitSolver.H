#ifndef _IMPLICIT_SOLVER_H_
#define _IMPLICIT_SOLVER_H_

#include "FieldSolver/ImplicitSolvers/Utils/WarpXSolverVec.H"

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

/*
 * Base class for implicit solvers. The base functions are those
 * needed by an implicit solver to be used through WarpX and those
 * needed to use the nonlinear solvers, such as Picard or Newton.
 */

class WarpX;
class ImplicitSolver
{
public:

    ImplicitSolver() {}

    virtual ~ImplicitSolver() = default;

    //
    // the following routines are called by WarpX
    //

    virtual void Define ( WarpX* const  a_WarpX ) = 0;
    
    virtual bool IsDefined () const = 0;
    
    virtual void PrintParameters () const = 0;

    virtual void Initialize () = 0;
    
    virtual void GetParticleSolverParams (int&  a_max_particle_iter,
                                          amrex::ParticleReal&  a_particle_tol ) = 0;
    
    virtual void OneStep ( const amrex::Real  a_time,
                           const amrex::Real  a_dt,
                           const int          a_step ) = 0;

    //
    // the following routines are called by the linear and nonlinear solvers
    //

    virtual void PreRHSOp ( const WarpXSolverVec&  a_E,
                            const amrex::Real      a_time,
                            const amrex::Real      a_dt,
                            const int              a_nl_iter,
                            const bool             a_from_jacobian ) = 0;

    virtual void ComputeRHS ( WarpXSolverVec&  a_Erhs,
                        const WarpXSolverVec&  a_E,
                        const amrex::Real      a_time,
                        const amrex::Real      a_dt ) = 0;

private:

};

#endif
