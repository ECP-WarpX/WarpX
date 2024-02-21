#ifndef _IMPLICIT_SOLVER_H_
#define _IMPLICIT_SOLVER_H_

#include "FieldSolver/ImplicitSolvers/Utils/WarpXSolverVec.H"

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

class WarpX;
class ImplicitSolver
{
    public:

    ImplicitSolver() {}

    virtual ~ImplicitSolver() = default;
    
    virtual void Define( WarpX* const  a_WarpX ) = 0;
    
    virtual bool isDefined() const = 0;
    
    virtual void PrintParameters() const = 0;

    virtual void Initialize() = 0;
    
    virtual void OneStep( const amrex::Real  a_time,
                          const amrex::Real  a_dt,
                          const int          a_step ) = 0;

    virtual void PreRHSOp( const WarpXSolverVec&  a_E,
                           const amrex::Real      a_time,
                           const amrex::Real      a_dt,
                           const int              a_nl_iter,
                           const bool             a_from_jacobian ) = 0;

    virtual void ComputeRHS( WarpXSolverVec&  a_Erhs,
                       const WarpXSolverVec&  a_E,
                       const amrex::Real      a_time,
                       const amrex::Real      a_dt ) = 0;

    virtual void PostUpdateState( const WarpXSolverVec&  a_E,
                                  const amrex::Real      a_time,
                                  const amrex::Real      a_dt ) = 0;

    private:

};

#endif
