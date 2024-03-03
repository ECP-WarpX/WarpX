/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef THETA_IMPLICIT_EM_H_
#define THETA_IMPLICIT_EM_H_

#include "FieldSolver/ImplicitSolvers/WarpXSolverVec.H"
#include "NonlinearSolvers/NonlinearSolverLibrary.H"

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "ImplicitSolver.H"

/** @file
 *  Electromagnetic theta-implicit time solver class. This is a fully implicit
 *  algorithm where both the fields and particles are treated implicitly.
 *
 *  The time stencil is
 *  Eg^{n+1} = Eg^{n} + c^2*dt*( curlBg^{n+theta} - mu0*Jg^{n+1/2} )
 *  Bg^{n+1} = Bg^{n} - dt*curlEg^{n+theta}
 *  xp^{n+1} = xp^n + dt*up^{n+1/2}/(gammap^n + gammap^{n+1})
 *  up^{n+1} = up^n + dt*qp/mp*(Ep^{n+theta} + up^{n+1/2}/gammap^{n+1/2} x Bp^{n+theta})
 *  where f^{n+theta} = (1.0-theta)*f^{n} + theta*f^{n+1} with 0.5 <= theta <= 1.0
 *
 *  The user-specified time-biasing parameter theta used for the fields on the RHS is bound
 *  between 0.5 and 1.0. The algorithm is exactly energy conserving for theta = 0.5.
 *  Signifcant damping of high-k modes will occur as theta approaches 1.0. The algorithm is
 *  numerially stable for any time step. I.e., the CFL condition for light waves does not
 *  have to be satisifed and the time step is not limited by the plasma period. However, how
 *  efficiently the algorithm can use large time steps depends strongly on the nonlinear solver.
 *  Furthermore, the time step should always be such that particles do not travel outside the
 *  ghost region of the box they live in, which is an MPI-related limitation. The time step
 *  is always limited by the need to resolve the appropriate physics.
 *
 *  See S. Markidis, G. Lapenta, "The energy conserving particle-in-cell method." JCP 230 (2011).
 *
 *  See G. Chen, L. Chacon, D.C. Barnes, "An energy- and charge-conserving, implicit,
 *  elctrostatic particle-in-cell algorithm." JCP 230 (2011).
 *
 *  See J.R. Angus, A. Link, A. Friedman, D. Ghosh, J. D. Johnson, "On numerical energy
 *  conservation for an implicit particle-in-cell method coupled with a binary Monte-Carlo
 *  algorithm for Coulomb collisions.", JCP 456 (2022).
 *
 *  See J.R. Angus, W. Farmer, A. Friedman, D. Ghosh, D. Grote, D. Larson, A. Link, "An
 *  implicit particle code with exact energy and charge consevation for electromagnetic studies
 *  of dense plasmas.", JCP 491 (2023).
 */

class ThetaImplicitEM : public ImplicitSolver
{
public:

    ThetaImplicitEM() = default;

    ~ThetaImplicitEM() override {
        // If we don't call the function below, then we get the
        // following error message after the simulation finishes:
        // malloc_consolidate(): unaligned fastbin chunk detected
        WarpXSolverVec::clearDotMask();
    }

    // Prohibit Move and Copy operations
    ThetaImplicitEM(const ThetaImplicitEM&) = delete;
    ThetaImplicitEM& operator=(const ThetaImplicitEM&) = delete;
    ThetaImplicitEM(ThetaImplicitEM&&) = delete;
    ThetaImplicitEM& operator=(ThetaImplicitEM&&) = delete;

    void Define ( WarpX*  a_WarpX ) override;

    [[nodiscard]] bool IsDefined () const override { return m_is_defined; }

    void PrintParameters () const override;

    void Initialize () override {;}

    void GetParticleSolverParams ( int&  a_max_particle_iterations,
                                   amrex::ParticleReal&  a_particle_tolerance ) override
    {
        a_max_particle_iterations = m_max_particle_iterations;
        a_particle_tolerance = m_particle_tolerance;
    }

    void OneStep ( amrex::Real  a_time,
                   amrex::Real  a_dt,
                   int          a_step ) override;

    void PreRHSOp ( const WarpXSolverVec&  a_E,
                          amrex::Real      a_time,
                          amrex::Real      a_dt,
                          int              a_nl_iter,
                          bool             a_from_jacobian ) override;

    void ComputeRHS ( WarpXSolverVec&  a_Erhs,
                const WarpXSolverVec&  a_E,
                      amrex::Real      a_time,
                      amrex::Real      a_dt ) override;

    [[nodiscard]] amrex::Real theta () const { return m_theta; }

private:

    bool m_is_defined = false;

    WarpX* m_WarpX;

    amrex::Real m_theta = 0.5;
    amrex::ParticleReal m_particle_tolerance = 1.0e-10;
    int m_max_particle_iterations = 20;

    NonlinearSolverType m_nlsolver_type = NonlinearSolverType::Picard;
    std::unique_ptr<NonlinearSolver<WarpXSolverVec,ThetaImplicitEM>> m_nlsolver;

    WarpXSolverVec m_E, m_Eold;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_Bold;

    void UpdateWarpXState ( const WarpXSolverVec&  a_E,
                            amrex::Real            a_time,
                            amrex::Real            a_dt );

};

#endif
