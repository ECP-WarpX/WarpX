
#include "ImplicitSolverEM.H"
#include "WarpX.H"

void ImplicitSolverEM::Define( WarpX* const  a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "ImplicitSolverEM object is already defined!");
    
    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;
    
    // initialize E vectors
    m_E.Define( m_WarpX->Efield_fp );
    m_Erhs  = m_E;
    m_Eold  = m_E;
    m_Esave = m_E;
    
    // initialize B vectors
    m_B.Define( m_WarpX->Bfield_fp );
    m_Brhs = m_B;
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) { 
        m_Bold  = m_B;
    }
    
    // parse implicit solver parameters
    amrex::ParmParse pp("algo.implicit");
    pp.query("verbose", m_verbose);
    pp.query("theta", m_theta);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_theta>=0.5 && m_theta<=1.0,
        "theta parameter for theta implicit time solver must be between 0.5 and 1.0");

    std::string nlsolver_type_str;
    pp.query("nonlinear_solver", nlsolver_type_str);
    if (nlsolver_type_str=="picard" || nlsolver_type_str=="Picard") {
        m_nlsolver_type = NonlinearSolverType::Picard;
    }
    else if (nlsolver_type_str=="newton" || nlsolver_type_str=="Newton") {
        m_nlsolver_type = NonlinearSolverType::Newton;
    }

    // parse picard solver parameters
    amrex::ParmParse pp_picard("algo.picard");
    pp_picard.query("require_convergence", m_require_convergence);
    pp_picard.query("relative_tolerance", m_rtol);
    pp_picard.query("max_iterations", m_max_iter);

    //if (m_nlsolver_type == _NLSOLVER_PICARD_) {

    //m_nlsolver = new PicardSolver<ODEVector<EMFields>, PICTimeIntegrator>;
    //m_nlsolver->define(m_U, this, m_func, m_theta);
    //m_nlsolver->setOutputIndent("  ");
  
    //Real rtol, atol;
    //int maxits;
    //m_nlsolver->getParams(rtol,atol,maxits);    
    //if( (m_nlsolver_type == _NLSOLVER_PICARD_) && (maxits % 2 == 0) ) {
    // cout << "WARNING!!! iter_max = " << maxits << endl;
    //  cout << "theta_implicit time advance with picard solver " << endl;
    //  cout << "may not work properly with even iter_max " << endl;
    //} 
  
    if (m_verbose) { PrintParams(); }
    m_is_defined = true;
}

void ImplicitSolverEM::PrintParams() const
{
    amrex::Print(0) << std::endl;
    amrex::Print(0) << "================== Implicit Solver ==================" << std::endl;
    amrex::Print(0)     << "time-bias parameter theta:  " << m_theta << std::endl;
    if (m_nlsolver_type==NonlinearSolverType::Picard) {
        amrex::Print(0) << "nonlinear solver type:      Picard" << std::endl;
    }
    else if (m_nlsolver_type==NonlinearSolverType::Newton) {
        amrex::Print(0) << "nonlinear solver type:      Newton" << std::endl;
    }
    amrex::Print(0)     << "picard require convergence: " << (m_require_convergence?"true":"false") << std::endl;
    amrex::Print(0)     << "picard relative tolerance:  " << m_rtol << std::endl;
    amrex::Print(0)     << "picard max iterations:      " << m_max_iter << std::endl;
    //cout << "solver_type = " << m_nlsolver_type << endl;
    //m_nlsolver->printParams();
    amrex::Print(0) << "=====================================================" << std::endl;
    amrex::Print(0) << std::endl;
}

void ImplicitSolverEM::Initialize( )
{

    // initialize E vectors
    m_E.Copy( m_WarpX->Efield_fp );
    m_Eold = m_E;
    
    // initialize B vectors
    m_B.Copy( m_WarpX->Bfield_fp );
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) { 
        m_Bold = m_B;
    }

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
    
    //const amrex::Real half_time = a_old_time + a_dt/2.0;

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save the fields at the start of the step
    //m_Eold = m_E;
    m_Eold.Copy(m_WarpX->Efield_fp);

    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        //m_Bold = m_B;
        m_Bold.Copy(m_WarpX->Bfield_fp);
    } else if (m_WarpX->evolve_scheme == EvolveScheme::SemiImplicit) {
        // Compute Bfield at time n+1/2
        m_WarpX->ComputeRHSB(a_dt, m_Brhs);
        m_B += m_Brhs;
        m_WarpX->UpdateMagneticField( m_B );
        //m_WarpX->Bfield_fp[0][0]->plus(*(m_Brhs.getVec()[0][0]), 0, 1, 0);
        //m_WarpX->Bfield_fp[0][1]->plus(*(m_Brhs.getVec()[0][1]), 0, 1, 0);
        //m_WarpX->Bfield_fp[0][2]->plus(*(m_Brhs.getVec()[0][2]), 0, 1, 0);
        //m_WarpX->ApplyMagneticFieldBCs( );
    }

    // Start the iterations
    amrex::Real norm, norm0;
    int iteration_count = 0;
    while (iteration_count < m_max_iter) {
        iteration_count++;

        // Advance the particle positions by 1/2 dt,
        // particle velocities by dt, then take average of old and new v,
        // deposit currents, giving J at n+1/2 used in ComputeRHSE below
        m_WarpX->PreRHSOpFromNonlinearIter( a_old_time, a_dt, iteration_count );
        
        // Save E at n+1/2 from the previous iteration to compute step norm
        //m_Esave = m_E;
        m_Esave.Copy(m_WarpX->Efield_fp);
        
        // Compute Efield at time n+1/2
        m_E = m_Eold;
        m_WarpX->ComputeRHSE(0.5_rt*a_dt, m_Erhs);
        m_E += m_Erhs;
        m_WarpX->UpdateElectricField( m_E );
        
        if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
            // Compute Bfield at time n+1/2
            m_B = m_Bold;
            m_WarpX->ComputeRHSB(0.5_rt*a_dt, m_Brhs);
            m_B += m_Brhs;
            m_WarpX->UpdateMagneticField( m_B );
        }

        // The B field update needs
        if (m_WarpX->num_mirrors>0){
            m_WarpX->applyMirrors(a_old_time);
            // E : guard cells are NOT up-to-date from the mirrors
            // B : guard cells are NOT up-to-date from the mirrors
        }
        
        // Compute the step norm for E and for B    
        m_Esave.getVec()[0][0]->minus(*m_WarpX->Efield_fp[0][0], 0, 1, 0);
        m_Esave.getVec()[0][1]->minus(*m_WarpX->Efield_fp[0][1], 0, 1, 0);
        m_Esave.getVec()[0][2]->minus(*m_WarpX->Efield_fp[0][2], 0, 1, 0);
        norm = m_Esave.norm();
        if (iteration_count==1) { 
            if (norm > 0.) { norm0 = norm; }
            else { norm0 = 1._rt; }
        }
            
        amrex::Real rnorm_E = norm/norm0;
        if (m_verbose || iteration_count == m_max_iter) {
            amrex::Print() << "Picard: iter = " << std::setw(3) << iteration_count <<  ", norm = " 
                           << std::scientific << std::setprecision(5) << norm << " (abs.), " 
                           << std::scientific << std::setprecision(5) << rnorm_E << " (rel.)" << "\n";
        }

        if (rnorm_E < m_rtol) {
            amrex::Print() << "Picard: exiting at iter = " << std::setw(3) << iteration_count 
                           << ". Satisified relative tolerance " << m_rtol << std::endl;
            break;
        }

    }

    if (m_rtol > 0. && iteration_count == m_max_iter) {
       std::stringstream convergenceMsg;
       convergenceMsg << "Picard solver failed to converge after " << iteration_count << 
                         " iterations. Relative norm is " << norm/norm0 << 
                         " and the relative tolerance is " << m_rtol;
       if (m_verbose) { amrex::Print() << convergenceMsg.str() << std::endl; }
       if (m_require_convergence) {
           WARPX_ABORT_WITH_MESSAGE(convergenceMsg.str());
       } else {
           ablastr::warn_manager::WMRecordWarning("PicardSolver", convergenceMsg.str());
       }
    }

    // update field boundary probes, and then advance fields to new time
    //m_fields->updateBoundaryProbes( a_dt );
    const amrex::Real new_time = a_old_time + a_dt;
    //m_U = (1.0/m_theta)*m_U + ((m_theta-1.0)/m_theta)*m_Uold;
    //m_fields->updatePhysicalState( m_U, new_time, e_and_b );
    
    // Advance particles to step n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance fields to step n+1
    m_WarpX->FinishImplicitFieldUpdate(m_WarpX->Efield_fp, m_Eold.getVec());
    m_WarpX->ApplyElectricFieldBCs( false );
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        m_WarpX->FinishImplicitFieldUpdate(m_WarpX->Bfield_fp, m_Bold.getVec());
        m_WarpX->ApplyMagneticFieldBCs( false );
    }
  
}

void ImplicitSolverEM::ComputeRHS( WarpXFieldVec&  a_Erhs_vec,
                             const WarpXFieldVec&  a_Efield_vec,
                             const amrex::Real     a_time,
                             const amrex::Real     a_dt )
{  
  
  // this function is called from the nonlinear solvers

}

void ImplicitSolverEM::UpdatePhysicalState( WarpXFieldVec&  a_field_vec,
                                      const amrex::Real     a_time )
{

}
