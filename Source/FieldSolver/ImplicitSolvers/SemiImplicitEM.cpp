/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SemiImplicitEM.H"
#include "WarpX.H"

using namespace warpx::fields;
using namespace amrex::literals;

void SemiImplicitEM::Define ( WarpX*  a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "SemiImplicitEM object is already defined!");

    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;

    // Define E and Eold vectors
    m_E.Define( m_WarpX->getMultiLevelField(FieldType::Efield_fp) );
    m_Eold.Define( m_WarpX->getMultiLevelField(FieldType::Efield_fp) );

    // Need to define the WarpXSolverVec owned dot_mask to do dot
    // product correctly for linear and nonlinear solvers
    const amrex::Vector<amrex::Geometry>& Geom = m_WarpX->Geom();
    m_E.SetDotMask(Geom);

    // Parse implicit solver parameters
    const amrex::ParmParse pp("implicit_evolve");
    parseNonlinearSolverParams( pp );

    // Define the nonlinear solver
    m_nlsolver->Define(m_E, this);
    m_is_defined = true;

}

void SemiImplicitEM::PrintParameters () const
{
    if (!m_WarpX->Verbose()) { return; }
    amrex::Print() << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "----------- SEMI IMPLICIT EM SOLVER PARAMETERS ------------" << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "max particle iterations:    " << m_max_particle_iterations << std::endl;
    amrex::Print() << "particle tolerance:         " << m_particle_tolerance << std::endl;
    if (m_nlsolver_type==NonlinearSolverType::Picard) {
        amrex::Print() << "Nonlinear solver type:      Picard" << std::endl;
    }
    else if (m_nlsolver_type==NonlinearSolverType::Newton) {
        amrex::Print() << "Nonlinear solver type:      Newton" << std::endl;
    }
    m_nlsolver->PrintParams();
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << std::endl;
}

void SemiImplicitEM::OneStep ( amrex::Real  a_time,
                               amrex::Real  a_dt,
                               int          a_step )
{
    amrex::ignore_unused(a_step);

    // Fields have E^{n}, B^{n-1/2}
    // Particles have p^{n} and x^{n}.

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save the fields at the start of the step
    m_Eold.Copy( m_WarpX->getMultiLevelField(FieldType::Efield_fp) );
    m_E.Copy(m_Eold); // initial guess for E

    // Compute Bfield at time n+1/2
    m_WarpX->EvolveB(a_dt, DtType::Full);
    m_WarpX->ApplyMagneticFieldBCs();

    const amrex::Real half_time = a_time + 0.5_rt*a_dt;

    // Solve nonlinear system for E at t_{n+1/2}
    // Particles will be advanced to t_{n+1/2}
    m_nlsolver->Solve( m_E, m_Eold, half_time, a_dt );

    // Update WarpX owned Efield_fp to t_{n+1/2}
    m_WarpX->SetElectricFieldAndApplyBCs( m_E );

    // Advance particles from time n+1/2 to time n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance E fields from time n+1/2 to time n+1
    // Eg^{n+1} = 2.0*E_g^{n+1/2} - E_g^n
    m_E.linComb( 2._rt, m_E, -1._rt, m_Eold );
    m_WarpX->SetElectricFieldAndApplyBCs( m_E );

}

void SemiImplicitEM::ComputeRHS ( WarpXSolverVec&  a_RHS,
                            const WarpXSolverVec&  a_E,
                                  amrex::Real      a_time,
                                  amrex::Real      a_dt,
                                  int              a_nl_iter,
                                  bool             a_from_jacobian )
{
    // update WarpX-owned Efield_fp using current state of E from
    // the nonlinear solver at time n+theta
    m_WarpX->SetElectricFieldAndApplyBCs( a_E );

    // Self consistently update particle positions and velocities using the
    // current state of the fields E and B. Deposit current density at time n+1/2.
    m_WarpX->ImplicitPreRHSOp( a_time, a_dt, a_nl_iter, a_from_jacobian );

    // RHS = cvac^2*0.5*dt*( curl(B^{n+1/2}) - mu0*J^{n+1/2} )
    m_WarpX->ImplicitComputeRHSE(0.5_rt*a_dt, a_RHS);
}
