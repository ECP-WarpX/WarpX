#ifndef _SEMI_IMPLICIT_EM_H_
#define _SEMI_IMPLICIT_EM_H_

#include "FieldSolver/ImplicitSolvers/Utils/WarpXSolverVec.H"
#include "NonlinearSolvers/NonlinearSolverLibrary.H"

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "ImplicitSolver.H"

/*
 *  semi-implicit electromagnetic time solver class.
 */

class SemiImplicitEM : public ImplicitSolver
{
public:

    SemiImplicitEM()
    {
        m_nlsolver_type = NonlinearSolverType::Picard;
        m_max_particle_iterations = 21;
        m_particle_tolerance = 1.0e-10;
        m_is_defined = false;
        m_verbose = false;
    }

    virtual ~SemiImplicitEM() {
        // If we don't call the function below, then we get the 
        // following error message after the simulation finishes:
        // malloc_consolidate(): unaligned fastbin chunk detected
        WarpXSolverVec::clearDotMask();
    }

    virtual void Define ( WarpX* const  a_WarpX );

    virtual bool IsDefined () const { return m_is_defined; }

    virtual void PrintParameters () const;

    virtual void Initialize () {;}

    virtual void GetParticleSolverParams (int&  a_max_particle_iterations,
                                          amrex::ParticleReal&  a_particle_tolerance )
    {
        a_max_particle_iterations = m_max_particle_iterations;
        a_particle_tolerance = m_particle_tolerance;
    }

    virtual void OneStep ( const amrex::Real  a_time,
                           const amrex::Real  a_dt,
                           const int          a_step );

    virtual void PreRHSOp ( const WarpXSolverVec&  a_E,
                            const amrex::Real      a_time,
                            const amrex::Real      a_dt,
                            const int              a_nl_iter,
                            const bool             a_from_jacobian );

    virtual void ComputeRHS ( WarpXSolverVec&  a_Erhs,
                        const WarpXSolverVec&  a_E,
                        const amrex::Real      a_time,
                        const amrex::Real      a_dt );

private:

    bool m_is_defined;
    bool m_verbose;

    WarpX* m_WarpX;

    amrex::ParticleReal m_particle_tolerance;
    int m_max_particle_iterations;

    NonlinearSolverType m_nlsolver_type;
    std::unique_ptr<NonlinearSolver<WarpXSolverVec,SemiImplicitEM>> m_nlsolver;

    WarpXSolverVec m_E, m_Eold;

    void UpdateWarpXState ( const WarpXSolverVec&  a_E,
                            const amrex::Real      a_time,
                            const amrex::Real      a_dt );

};

#endif
