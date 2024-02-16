
#include "ImplicitSolverEM.H"
#include "WarpX.H"

void ImplicitSolverEM::Define( WarpX* const  a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "ImplicitSolverEM object is already defined!");
    
    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;
    
    // Define E vectors
    m_E.Define( m_WarpX->Efield_fp );
    m_Eold  = m_E;
    
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) { 
        // Define Bold vector
        const int lev = 0;
        m_Bold.resize(1); // size is number of levels
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& Bmf = *(m_WarpX->Bfield_fp[lev][n]);
            m_Bold[lev][n] = std::make_unique<amrex::MultiFab>( Bmf.boxArray(), 
                                                                Bmf.DistributionMap(),
                                                                Bmf.nComp(), 
                                                                Bmf.nGrowVect() );
        }
    }
    
    // Parse implicit solver parameters
    amrex::ParmParse pp("implicit_evolve");
    pp.query("verbose", m_verbose);
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) { 
        pp.query("theta", m_theta);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_theta>=0.5 && m_theta<=1.0,
            "theta parameter for theta implicit time solver must be between 0.5 and 1.0");
    }

    std::string nlsolver_type_str;
    pp.query("nonlinear_solver", nlsolver_type_str);
    if (nlsolver_type_str=="picard" || nlsolver_type_str=="Picard") {
        m_nlsolver_type = NonlinearSolverType::Picard;
    }
    else if (nlsolver_type_str=="newton" || nlsolver_type_str=="Newton") {
        m_nlsolver_type = NonlinearSolverType::Newton;
    }

    // Define the nonlinear solver
    if (m_nlsolver_type == NonlinearSolverType::Picard) {
        m_nlsolver = std::make_unique<PicardSolver<WarpXSolverVec,ImplicitSolverEM>>();
        m_nlsolver->Define(m_E, this);
    }

    //if (m_verbose) { PrintParameters(); }
    m_is_defined = true;
}

void ImplicitSolverEM::PrintParameters() const
{
    if (!m_verbose) { return; }
    amrex::Print() << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "-------------- IMPLICIT EM SOLVER PARAMETERS --------------" << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print()     << "Time-bias parameter theta:  " << m_theta << std::endl;
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

void ImplicitSolverEM::Initialize( )
{

    // initialize E vectors
    m_E.Copy( m_WarpX->Efield_fp );
    m_Eold = m_E;
    
}

void ImplicitSolverEM::OneStep( const amrex::Real  a_old_time,
                                const amrex::Real  a_dt,
                                const int          a_step )
{  
    using namespace amrex::literals;
    amrex::ignore_unused(a_step);

    // We have E^{n}.
    // Particles have p^{n} and x^{n}.
    // With full implicit, B^{n}
    // With semi-implicit, B^{n-1/2}

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save the fields at the start of the step
    m_Eold.Copy(m_WarpX->Efield_fp);
    m_E = m_Eold; // initial guess for E

    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        int lev = 0;
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& Bfp  = *(m_WarpX->Bfield_fp[lev][n]);
            amrex::MultiFab& Bold = *m_Bold[lev][n];
            amrex::MultiFab::Copy(Bold, Bfp, 0, 0, 1, Bold.nGrowVect());
        }
    } else if (m_WarpX->evolve_scheme == EvolveScheme::SemiImplicit) {
        // Compute Bfield at time n+1/2
        m_WarpX->EvolveB(a_dt, DtType::Full);
        m_WarpX->ApplyMagneticFieldBCs( true );
    }

    // Solve nonlinear system for E at t_{n+theta}
    // B will also be advanced to t_{n+theta}
    // Particles will be advanced to t_{n+1/2}
    m_nlsolver->Solve( m_E, m_Eold, a_old_time, a_dt );

    // Update field boundary probes prior to updating fields to t_{n+1}
    //m_fields->updateBoundaryProbes( a_dt );
    
    // Advance particles to step n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance fields to step n+1
    m_WarpX->FinishImplicitFieldUpdate(m_E.getVec(), m_Eold.getVec(), m_theta);
    m_WarpX->UpdateElectricField( m_E, false ); // JRA not sure about false here. is what DG had.
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        m_WarpX->FinishImplicitFieldUpdate(m_WarpX->Bfield_fp, m_Bold, m_theta);
        m_WarpX->ApplyMagneticFieldBCs( false );
    }
  
}

void ImplicitSolverEM::PreRHSOp( const WarpXSolverVec&  a_E,
                                 const amrex::Real      a_time,
                                 const amrex::Real      a_dt,
                                 const int              a_nl_iter )
{  
    amrex::ignore_unused(a_E);
    
    if (m_nlsolver_type!=NonlinearSolverType::Picard) {
        PostUpdateState( a_E, a_time, a_dt );
    }

    // Advance the particle positions by 1/2 dt,
    // particle velocities by dt, then take average of old and new v,
    // deposit currents, giving J at n+1/2 used in ComputeRHSE below
    m_WarpX->PreRHSOpFromNonlinearIter( a_time, a_dt, a_nl_iter );

}

void ImplicitSolverEM::ComputeRHS( WarpXSolverVec&  a_Erhs,
                             const WarpXSolverVec&  a_E,
                             const amrex::Real      a_time,
                             const amrex::Real      a_dt )
{  
    amrex::ignore_unused(a_E, a_time);
    m_WarpX->ComputeRHSE(m_theta*a_dt, a_Erhs);

}

void ImplicitSolverEM::PostUpdateState( const WarpXSolverVec&  a_E,
                                        const amrex::Real      a_time,
                                        const amrex::Real      a_dt )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_time);

    // Update Efield_fp owned by WarpX
    m_WarpX->UpdateElectricField( a_E, true );

    // Update Bfield_fp owned by WarpX
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        // Compute Bfield at time n+theta
        const int lev = 0;
        for (int n=0; n<3; n++) {
            amrex::MultiFab& Bfp  = *(m_WarpX->Bfield_fp[lev][n]);
            const amrex::MultiFab& Bold = *m_Bold[lev][n];
            amrex::MultiFab::Copy(Bfp, Bold, 0, 0, 1, Bold.nGrowVect());
        }
        m_WarpX->EvolveB(m_theta*a_dt, DtType::Full);
        m_WarpX->ApplyMagneticFieldBCs( true );
    }

    // The B field update needs. Talk to DG about this. Only needed when B updates?
    if (m_WarpX->num_mirrors>0){
        m_WarpX->applyMirrors(a_time);
        // E : guard cells are NOT up-to-date from the mirrors
        // B : guard cells are NOT up-to-date from the mirrors
    }

}
