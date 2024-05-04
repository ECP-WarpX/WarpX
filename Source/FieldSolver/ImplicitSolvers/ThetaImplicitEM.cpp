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
    m_E.Define( m_WarpX->getMultiLevelField(FieldType::Efield_fp) );
    m_Eold.Define( m_WarpX->getMultiLevelField(FieldType::Efield_fp) );

    // Need to define the WarpXSolverVec owned dot_mask to do dot
    // product correctly for linear and nonlinear solvers
    const amrex::Vector<amrex::Geometry>& Geom = m_WarpX->Geom();
    m_E.SetDotMask(Geom);

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
    amrex::Print() << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "----------- THETA IMPLICIT EM SOLVER PARAMETERS -----------" << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "Time-bias parameter theta:  " << m_theta << std::endl;
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

void ThetaImplicitEM::OneStep ( const amrex::Real  a_time,
                                const amrex::Real  a_dt,
                                const int          a_step )
{
    amrex::ignore_unused(a_step);

    // Fields have E^{n} and B^{n}
    // Particles have p^{n} and x^{n}.

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save the fields at the start of the step
    m_Eold.Copy( m_WarpX->getMultiLevelField(FieldType::Efield_fp) );
    m_E = m_Eold; // initial guess for E

    const int num_levels = m_Bold.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& Bfp = m_WarpX->getField(FieldType::Bfield_fp,lev,n);
            amrex::MultiFab& Bold = *m_Bold[lev][n];
            amrex::MultiFab::Copy(Bold, Bfp, 0, 0, 1, Bold.nGrowVect());
        }
    }

    const amrex::Real theta_time = a_time + m_theta*a_dt;

    // Solve nonlinear system for E at t_{n+theta}
    // Particles will be advanced to t_{n+1/2}
    m_nlsolver->Solve( m_E, m_Eold, theta_time, a_dt );

    // Update WarpX owned Efield_fp and Bfield_fp to t_{n+theta}
    UpdateWarpXFields( m_E, theta_time, a_dt );

    // Advance particles from time n+1/2 to time n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance E and B fields from time n+theta to time n+1
    const amrex::Real new_time = a_time + a_dt;
    FinishFieldUpdate( new_time );

}

void ThetaImplicitEM::PreRHSOp ( const WarpXSolverVec&  a_E,
                                 const amrex::Real      a_time,
                                 const amrex::Real      a_dt,
                                 const int              a_nl_iter,
                                 const bool             a_from_jacobian )
{

    // update derived variable B and then update WarpX owned Efield_fp and Bfield_fp
    UpdateWarpXFields( a_E, a_time, a_dt );

    // Advance the particle positions by 1/2 dt,
    // particle velocities by dt, then take average of old and new v,
    // deposit currents, giving J at n+1/2 used in ImplicitComputeRHSE below
    m_WarpX->ImplicitPreRHSOp( a_time, a_dt, a_nl_iter, a_from_jacobian );

}

void ThetaImplicitEM::ComputeRHS ( WarpXSolverVec&  a_Erhs,
                             const WarpXSolverVec&  a_E,
                                   amrex::Real      a_time,
                                   amrex::Real      a_dt )
{
    amrex::ignore_unused(a_E, a_time);
    m_WarpX->ImplicitComputeRHSE(m_theta*a_dt, a_Erhs);
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

    // Eg^{n+1} = (1/theta)*E_g^{n+theta} + (1-1/theta)*E_g^n
    // Bg^{n+1} = (1/theta)*B_g^{n+theta} + (1-1/theta)*B_g^n

    const amrex::Real c0 = 1._rt/m_theta;
    const amrex::Real c1 = 1._rt - c0;
    m_E.linComb( c0, m_E, c1, m_Eold );
    m_WarpX->SetElectricFieldAndApplyBCs( m_E );
    m_WarpX->FinishMagneticFieldAndApplyBCs( m_Bold, m_theta );

}
