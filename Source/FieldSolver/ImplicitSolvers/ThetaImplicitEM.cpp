/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ThetaImplicitEM.H"
#include "WarpX.H"

void ThetaImplicitEM::Define ( WarpX* const  a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "ThetaImplicitEM object is already defined!");

    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;

    // Define E vectors
    m_E.Define( m_WarpX->getEfield_fp_vec() );
    m_Eold.Define( m_WarpX->getEfield_fp_vec() );

    // Need to define the WarpXSolverVec owned dot_mask to do dot
    // product correctly for linear and nonlinear solvers
    const amrex::Vector<amrex::Geometry>& Geom = m_WarpX->Geom();
    m_E.SetDotMask(Geom);

    // Define Bold vector
    const int lev = 0;
    m_Bold.resize(1); // size is number of levels
    for (int n=0; n<3; n++) {
        const amrex::MultiFab& Bfp = m_WarpX->getBfield_fp(lev,n);
        m_Bold[lev][n] = std::make_unique<amrex::MultiFab>( Bfp.boxArray(),
                                                            Bfp.DistributionMap(),
                                                            Bfp.nComp(),
                                                            Bfp.nGrowVect() );
    }

    // Parse implicit solver parameters
    amrex::ParmParse pp("implicit_evolve");
    pp.query("verbose", m_verbose);
    pp.query("theta", m_theta);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_theta>=0.5 && m_theta<=1.0,
        "theta parameter for theta implicit time solver must be between 0.5 and 1.0");

    std::string nlsolver_type_str;
    pp.query("nonlinear_solver", nlsolver_type_str);
    if (nlsolver_type_str=="picard") {
        m_nlsolver_type = NonlinearSolverType::Picard;
        m_max_particle_iterations = 1;
        m_particle_tolerance = 0.0;
    }
    else if (nlsolver_type_str=="newton") {
        m_nlsolver_type = NonlinearSolverType::Newton;
        pp.query("max_particle_iterations", m_max_particle_iterations);
        pp.query("particle_tolerance", m_particle_tolerance);
    }
    else {
        WARPX_ABORT_WITH_MESSAGE(
            "invalid nonlinear_solver specified. Valid options are picard and newton.");
    }

    // Define the nonlinear solver
    if (m_nlsolver_type == NonlinearSolverType::Picard) {
        m_nlsolver = std::make_unique<PicardSolver<WarpXSolverVec,ThetaImplicitEM>>();
        m_nlsolver->Define(m_E, this);
    }
    else if (m_nlsolver_type == NonlinearSolverType::Newton) {
        m_nlsolver = std::make_unique<NewtonSolver<WarpXSolverVec,ThetaImplicitEM>>();
        m_nlsolver->Define(m_E, this);
    }

    m_is_defined = true;
}

void ThetaImplicitEM::PrintParameters () const
{
    if (!m_verbose) { return; }
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

void ThetaImplicitEM::OneStep ( const amrex::Real  a_old_time,
                                const amrex::Real  a_dt,
                                const int          a_step )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_step);

    // Fields have E^{n} and B^{n}
    // Particles have p^{n} and x^{n}.

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save the fields at the start of the step
    m_Eold.Copy( m_WarpX->getEfield_fp_vec() );
    m_E = m_Eold; // initial guess for E

    const int lev = 0;
    for (int n=0; n<3; n++) {
        const amrex::MultiFab& Bfp = m_WarpX->getBfield_fp(lev,n);
        amrex::MultiFab& Bold = *m_Bold[lev][n];
        amrex::MultiFab::Copy(Bold, Bfp, 0, 0, 1, Bold.nGrowVect());
    }

    // Solve nonlinear system for E at t_{n+theta}
    // Particles will be advanced to t_{n+1/2}
    m_nlsolver->Solve( m_E, m_Eold, a_old_time, a_dt );

    // update WarpX owned Efield_fp and Bfield_fp to t_{n+theta}
    UpdateWarpXState( m_E, a_old_time, a_dt );

    // Update field boundary probes prior to updating fields to t_{n+1}
    //m_fields->updateBoundaryProbes( a_dt );

    // Advance particles to step n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance fields to step n+1
    m_WarpX->FinishImplicitField(m_E.getVec(), m_Eold.getVec(), m_theta);
    m_WarpX->UpdateElectricField( m_E );
    m_WarpX->FinishMagneticField( m_Bold, m_theta );

}

void ThetaImplicitEM::PreRHSOp ( const WarpXSolverVec&  a_E,
                                 const amrex::Real      a_time,
                                 const amrex::Real      a_dt,
                                 const int              a_nl_iter,
                                 const bool             a_from_jacobian )
{
    amrex::ignore_unused(a_E);

    // update derived variable B and then update WarpX owned Efield_fp and Bfield_fp
    UpdateWarpXState( a_E, a_time, a_dt );

    // Advance the particle positions by 1/2 dt,
    // particle velocities by dt, then take average of old and new v,
    // deposit currents, giving J at n+1/2 used in ComputeRHSE below
    m_WarpX->PreRHSOp( a_time, a_dt, a_nl_iter, a_from_jacobian );

}

void ThetaImplicitEM::ComputeRHS ( WarpXSolverVec&  a_Erhs,
                             const WarpXSolverVec&  a_E,
                             const amrex::Real      a_time,
                             const amrex::Real      a_dt )
{
    amrex::ignore_unused(a_E, a_time);
    m_WarpX->ComputeRHSE(m_theta*a_dt, a_Erhs);
}

void ThetaImplicitEM::UpdateWarpXState ( const WarpXSolverVec&  a_E,
                                         const amrex::Real      a_time,
                                         const amrex::Real      a_dt )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_time);

    // Update Efield_fp owned by WarpX
    m_WarpX->UpdateElectricField( a_E );

    // Update Bfield owned by WarpX
    m_WarpX->UpdateMagneticField( m_Bold, m_theta*a_dt );

    // The B field update needs. Talk to DG about this. Only needed when B updates?
    if (m_WarpX->num_mirrors>0){
        m_WarpX->applyMirrors(a_time);
        // E : guard cells are NOT up-to-date from the mirrors
        // B : guard cells are NOT up-to-date from the mirrors
    }

}
