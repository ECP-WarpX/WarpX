/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef SEMI_IMPLICIT_EM_H_
#define SEMI_IMPLICIT_EM_H_

#include "FieldSolver/ImplicitSolvers/WarpXSolverVec.H"
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

    SemiImplicitEM() = default;

    virtual ~SemiImplicitEM() {
        // If we don't call the function below, then we get the
        // following error message after the simulation finishes:
        // malloc_consolidate(): unaligned fastbin chunk detected
        WarpXSolverVec::clearDotMask();
    }

    virtual void Define ( WarpX*  a_WarpX );

    virtual bool IsDefined () const { return m_is_defined; }

    virtual void PrintParameters () const;

    virtual void Initialize () {;}

    virtual void GetParticleSolverParams (int&  a_max_particle_iterations,
                                          amrex::ParticleReal&  a_particle_tolerance )
    {
        a_max_particle_iterations = m_max_particle_iterations;
        a_particle_tolerance = m_particle_tolerance;
    }

    virtual void OneStep ( amrex::Real  a_time,
                           amrex::Real  a_dt,
                           int          a_step );

    virtual void PreRHSOp ( const WarpXSolverVec&  a_E,
                            amrex::Real            a_time,
                            amrex::Real            a_dt,
                            int                    a_nl_iter,
                            bool                   a_from_jacobian );

    virtual void ComputeRHS ( WarpXSolverVec&  a_Erhs,
                        const WarpXSolverVec&  a_E,
                              amrex::Real      a_time,
                              amrex::Real      a_dt );

private:

    bool m_is_defined = false;

    WarpX* m_WarpX;

    amrex::ParticleReal m_particle_tolerance = 1.0e-10;
    int m_max_particle_iterations = 20;

    NonlinearSolverType m_nlsolver_type = NonlinearSolverType::Picard;
    std::unique_ptr<NonlinearSolver<WarpXSolverVec,SemiImplicitEM>> m_nlsolver;

    WarpXSolverVec m_E, m_Eold;

    void UpdateWarpXState ( const WarpXSolverVec&  a_E,
                            amrex::Real            a_time,
                            amrex::Real            a_dt );

};

#endif
