
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
    m_Erhs  = m_E;
    m_Eold  = m_E;
    m_Esave = m_E;
    
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
    
    // parse implicit solver parameters
    amrex::ParmParse pp("algo.implicit");
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

    // parse picard solver parameters
    amrex::ParmParse pp_picard("algo.picard");
    pp_picard.query("require_convergence", m_require_convergence);
    pp_picard.query("relative_tolerance", m_rtol);
    pp_picard.query("max_iterations", m_max_iter);

    //if (m_nlsolver_type == NonlinearSolverType::Picard) {
    //    m_nlsolver = new PicardSolver<ODEVector<EMFields>, PICTimeIntegrator>;
    //    m_nlsolver->define(m_U, this, m_func, m_theta);
    //}

    if (m_verbose) { PrintParams(); }
    m_is_defined = true;
}

void ImplicitSolverEM::PrintParams() const
{
    amrex::Print() << std::endl;
    amrex::Print() << "================== Implicit Solver ==================" << std::endl;
    amrex::Print()     << "time-bias parameter theta:  " << m_theta << std::endl;
    if (m_nlsolver_type==NonlinearSolverType::Picard) {
        amrex::Print() << "nonlinear solver type:      Picard" << std::endl;
    }
    else if (m_nlsolver_type==NonlinearSolverType::Newton) {
        amrex::Print() << "nonlinear solver type:      Newton" << std::endl;
    }
    amrex::Print()     << "picard require convergence: " << (m_require_convergence?"true":"false") << std::endl;
    amrex::Print()     << "picard relative tolerance:  " << m_rtol << std::endl;
    amrex::Print()     << "picard max iterations:      " << m_max_iter << std::endl;
    //cout << "solver_type = " << m_nlsolver_type << endl;
    //m_nlsolver->printParams();
    amrex::Print() << "=====================================================" << std::endl;
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
            amrex::MultiFab::Copy(Bold, Bfp, 0, 0, 1, m_Bold[lev][n]->nGrowVect());
        }
    } else if (m_WarpX->evolve_scheme == EvolveScheme::SemiImplicit) {
        // Compute Bfield at time n+1/2
        m_WarpX->EvolveB(a_dt, DtType::Full);
        m_WarpX->ApplyMagneticFieldBCs( true );
    }

    // Start the iterations
    amrex::Real norm, norm0;
    int iteration_count = 0;
    while (iteration_count < m_max_iter) {

        PreRHSOp( m_E, a_old_time, a_dt, iteration_count );
        
        // Save E from the previous iteration to compute step norm
        m_Esave = m_E;
        
        // Compute Efield at time n+theta
        m_WarpX->ComputeRHSE(m_theta*a_dt, m_Erhs);
        m_E = m_Eold + m_Erhs;
        m_WarpX->UpdateElectricField( m_E, true );
        m_E.Copy(m_WarpX->Efield_fp);

        PostUpdateState( m_E, a_old_time, a_dt );

        // Compute the step norm of E
        m_Esave -= m_E;
        norm = m_Esave.norm();
        if (iteration_count==0) { 
            if (norm > 0.) { norm0 = norm; }
            else { norm0 = 1._rt; }
        }
        iteration_count++;
            
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

void ImplicitSolverEM::PreRHSOp( const WarpXFieldVec&  a_Efield_vec,
                                 const amrex::Real     a_time,
                                 const amrex::Real     a_dt,
                                 const int             a_nl_iter )
{  
    amrex::ignore_unused(a_Efield_vec);
    
    if (m_nlsolver_type!=NonlinearSolverType::Picard) {
        m_WarpX->UpdateElectricField( a_Efield_vec, true );
        PostUpdateState( a_Efield_vec, a_time, a_dt );
    }

    // Advance the particle positions by 1/2 dt,
    // particle velocities by dt, then take average of old and new v,
    // deposit currents, giving J at n+1/2 used in ComputeRHSE below
    m_WarpX->PreRHSOpFromNonlinearIter( a_time, a_dt, a_nl_iter );

}

void ImplicitSolverEM::ComputeRHS( WarpXFieldVec&  a_Erhs_vec,
                             const WarpXFieldVec&  a_Efield_vec,
                             const amrex::Real     a_time,
                             const amrex::Real     a_dt )
{  
  
  // this function is called from the nonlinear solvers

}

void ImplicitSolverEM::UpdateState( WarpXFieldVec&  a_E,
                              const amrex::Real     a_time )
{

}

void ImplicitSolverEM::PostUpdateState( const WarpXFieldVec&  a_field_vec,
                                        const amrex::Real     a_time,
                                        const amrex::Real     a_dt )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_field_vec, a_time);

    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        // Compute Bfield at time n+theta
        const int lev = 0;
        for (int n=0; n<3; n++) {
            amrex::MultiFab& Bfp  = *(m_WarpX->Bfield_fp[lev][n]);
            const amrex::MultiFab& Bold = *m_Bold[lev][n];
            amrex::MultiFab::Copy(Bfp, Bold, 0, 0, 1, m_Bold[lev][n]->nGrowVect());
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
