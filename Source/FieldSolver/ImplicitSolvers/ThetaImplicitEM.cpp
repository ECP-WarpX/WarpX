/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldSolver/Fields.H"
#include "ThetaImplicitEM.H"
#include "WarpX.H"

using namespace warpx::fields;
using namespace amrex::literals;

void ThetaImplicitEM::Define ( WarpX* const  a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "ThetaImplicitEM object is already defined!");

    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;

    // Define E and Eold vectors
    m_E.Define( m_WarpX, FieldType::Efield_fp );
    m_Eold.Define( m_E );

    // Define Bold MultiFab
    const int num_levels = 1;
    m_Bold.resize(num_levels); // size is number of levels
    for (int lev = 0; lev < num_levels; ++lev) {
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& Bfp = m_WarpX->getField( FieldType::Bfield_fp,lev,n);
            m_Bold[lev][n] = std::make_unique<amrex::MultiFab>( Bfp.boxArray(),
                                                                Bfp.DistributionMap(),
                                                                Bfp.nComp(),
                                                                Bfp.nGrowVect() );
        }
    }

    // Parse theta-implicit solver specific parameters
    const amrex::ParmParse pp("implicit_evolve");
    pp.query("theta", m_theta);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_theta>=0.5 && m_theta<=1.0,
        "theta parameter for theta implicit time solver must be between 0.5 and 1.0");

    // Parse nonlinear solver parameters
    parseNonlinearSolverParams( pp );

    // Define the nonlinear solver
    m_nlsolver->Define(m_E, this);
    m_is_defined = true;

}

void ThetaImplicitEM::PrintParameters () const
{
    if (!m_WarpX->Verbose()) { return; }
    amrex::Print() << "\n";
    amrex::Print() << "-----------------------------------------------------------\n";
    amrex::Print() << "----------- THETA IMPLICIT EM SOLVER PARAMETERS -----------\n";
    amrex::Print() << "-----------------------------------------------------------\n";
    amrex::Print() << "Time-bias parameter theta:  " << m_theta << "\n";
    amrex::Print() << "max particle iterations:    " << m_max_particle_iterations << "\n";
    amrex::Print() << "particle tolerance:         " << m_particle_tolerance << "\n";
    if (m_nlsolver_type==NonlinearSolverType::Picard) {
        amrex::Print() << "Nonlinear solver type:      Picard\n";
    }
    else if (m_nlsolver_type==NonlinearSolverType::Newton) {
        amrex::Print() << "Nonlinear solver type:      Newton\n";
    }
    m_nlsolver->PrintParams();
    amrex::Print() << "-----------------------------------------------------------\n\n";
}

void ThetaImplicitEM::OneStep ( const amrex::Real  a_time,
                                const amrex::Real  a_dt,
                                const int          a_step )
{
    amrex::ignore_unused(a_step);

    // Fields have Eg^{n} and Bg^{n}
    // Particles have up^{n} and xp^{n}.

    // Save up and xp at the start of the time step
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save Eg at the start of the time step
    m_Eold.Copy( FieldType::Efield_fp );

    const int num_levels = static_cast<int>(m_Bold.size());
    for (int lev = 0; lev < num_levels; ++lev) {
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& Bfp = m_WarpX->getField(FieldType::Bfield_fp,lev,n);
            amrex::MultiFab& Bold = *m_Bold[lev][n];
            amrex::MultiFab::Copy(Bold, Bfp, 0, 0, 1, Bold.nGrowVect());
        }
    }

    const amrex::Real theta_time = a_time + m_theta*a_dt;

    // Solve nonlinear system for Eg at t_{n+theta}
    // Particles will be advanced to t_{n+1/2}
    m_E.Copy(m_Eold); // initial guess for Eg^{n+theta}
    m_nlsolver->Solve( m_E, m_Eold, theta_time, a_dt );

    // Update WarpX owned Efield_fp and Bfield_fp to t_{n+theta}
    UpdateWarpXFields( m_E, theta_time, a_dt );

    // Advance particles from time n+1/2 to time n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance Eg and Bg from time n+theta to time n+1
    const amrex::Real new_time = a_time + a_dt;
    FinishFieldUpdate( new_time );

}

void ThetaImplicitEM::ComputeRHS ( WarpXSolverVec&  a_RHS,
                             const WarpXSolverVec&  a_E,
                                   amrex::Real      a_time,
                                   amrex::Real      a_dt,
                                   int              a_nl_iter,
                                   bool             a_from_jacobian )
{
    // Update WarpX-owned Efield_fp and Bfield_fp using current state of
    // Eg from the nonlinear solver at time n+theta
    UpdateWarpXFields( a_E, a_time, a_dt );

    // Update particle positions and velocities using the current state
    // of Eg and Bg. Deposit current density at time n+1/2
    m_WarpX->ImplicitPreRHSOp( a_time, a_dt, a_nl_iter, a_from_jacobian );

    // RHS = cvac^2*m_theta*dt*( curl(Bg^{n+theta}) - mu0*Jg^{n+1/2} )
    m_WarpX->ImplicitComputeRHSE(m_theta*a_dt, a_RHS);
}

void ThetaImplicitEM::UpdateWarpXFields ( const WarpXSolverVec&  a_E,
                                          amrex::Real            a_time,
                                          amrex::Real            a_dt )
{
    amrex::ignore_unused(a_time);

    // Update Efield_fp owned by WarpX
    m_WarpX->SetElectricFieldAndApplyBCs( a_E );

    // Update Bfield_fp owned by WarpX
    m_WarpX->UpdateMagneticFieldAndApplyBCs( m_Bold, m_theta*a_dt );

}

void ThetaImplicitEM::FinishFieldUpdate ( amrex::Real  a_new_time )
{
    amrex::ignore_unused(a_new_time);

    // Eg^{n+1} = (1/theta)*Eg^{n+theta} + (1-1/theta)*Eg^n
    // Bg^{n+1} = (1/theta)*Bg^{n+theta} + (1-1/theta)*Bg^n

    const amrex::Real c0 = 1._rt/m_theta;
    const amrex::Real c1 = 1._rt - c0;
    m_E.linComb( c0, m_E, c1, m_Eold );
    m_WarpX->SetElectricFieldAndApplyBCs( m_E );
    m_WarpX->FinishMagneticFieldAndApplyBCs( m_Bold, m_theta );

}
