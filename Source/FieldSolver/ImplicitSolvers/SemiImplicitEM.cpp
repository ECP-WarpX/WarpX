
#include "SemiImplicitEM.H"
#include "WarpX.H"

void SemiImplicitEM::Define ( WarpX* const  a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "SemiImplicitEM object is already defined!");
    
    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;
 
    // Define E vectors
    m_E.Define( m_WarpX->getEfield_fp_vec() );
    m_Eold.Define( m_WarpX->getEfield_fp_vec() );
    
    // Need to define the WarpXSolverVec owned dot_mask to do dot
    // product correctly for linear and nonlinear solvers
    const amrex::Vector<amrex::Geometry>& Geom = m_WarpX->Geom();
    m_E.SetDotMask(Geom);
    
    // Parse implicit solver parameters
    amrex::ParmParse pp("implicit_evolve");
    pp.query("verbose", m_verbose);

    std::string nlsolver_type_str;
    pp.query("nonlinear_solver", nlsolver_type_str);
    if (nlsolver_type_str=="picard") {
        m_nlsolver_type = NonlinearSolverType::Picard;
        m_max_particle_iterations = 1;
        m_particle_tolerance = 0.0;
        m_nlsolver = std::make_unique<PicardSolver<WarpXSolverVec,SemiImplicitEM>>();
    }
    else if (nlsolver_type_str=="newton") {
        m_nlsolver_type = NonlinearSolverType::Newton;
        pp.query("max_particle_iterations", m_max_particle_iterations);
        pp.query("particle_tolerance", m_particle_tolerance);
        m_nlsolver = std::make_unique<NewtonSolver<WarpXSolverVec,SemiImplicitEM>>();
    }
    else {
        WARPX_ABORT_WITH_MESSAGE(
            "invalid nonlinear_solver specified. Valid options are picard and newton.");
    }

    // Define the nonlinear solver
    m_nlsolver->Define(m_E, this);

    m_is_defined = true;
}

void SemiImplicitEM::PrintParameters () const
{
    if (!m_verbose) { return; }
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

void SemiImplicitEM::OneStep ( const amrex::Real  a_old_time,
                               const amrex::Real  a_dt,
                               const int          a_step )
{  
    using namespace amrex::literals;
    amrex::ignore_unused(a_step);

    // Fields have E^{n}, B^{n-1/2]
    // Particles have p^{n} and x^{n}.

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save the fields at the start of the step
    m_Eold.Copy( m_WarpX->getEfield_fp_vec() );
    m_E = m_Eold; // initial guess for E

    // Compute Bfield at time n+1/2
    m_WarpX->EvolveB(a_dt, DtType::Full);
    m_WarpX->ApplyMagneticFieldBCs( true );

    // Solve nonlinear system for E at t_{n+1/2}
    // Particles will be advanced to t_{n+1/2}
    m_nlsolver->Solve( m_E, m_Eold, a_old_time, a_dt );
    
    // update WarpX owned Efield_fp and Bfield_fp to t_{n+1/2}
    UpdateWarpXState( m_E, a_old_time, a_dt );

    // Update field boundary probes prior to updating fields to t_{n+1}
    //m_fields->updateBoundaryProbes( a_dt );
    
    // Advance particles to step n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance fields to step n+1
    m_WarpX->FinishImplicitField(m_E.getVec(), m_Eold.getVec(), 0.5);
    m_WarpX->UpdateElectricField( m_E, false ); // JRA not sure about false here. is what DG had.
  
}

void SemiImplicitEM::PreRHSOp ( const WarpXSolverVec&  a_E,
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

void SemiImplicitEM::ComputeRHS ( WarpXSolverVec&  a_Erhs,
                            const WarpXSolverVec&  a_E,
                            const amrex::Real      a_time,
                            const amrex::Real      a_dt )
{  
    amrex::ignore_unused(a_E, a_time);
    m_WarpX->ComputeRHSE(0.5*a_dt, a_Erhs);
}

void SemiImplicitEM::UpdateWarpXState ( const WarpXSolverVec&  a_E,
                                        const amrex::Real      a_time,
                                        const amrex::Real      a_dt )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_time,a_dt);

    // Update Efield_fp owned by WarpX
    m_WarpX->UpdateElectricField( a_E, true );

    // The B field update needs. Talk to DG about this. Only needed when B updates?
    if (m_WarpX->num_mirrors>0){
        m_WarpX->applyMirrors(a_time);
        // E : guard cells are NOT up-to-date from the mirrors
        // B : guard cells are NOT up-to-date from the mirrors
    }

}
