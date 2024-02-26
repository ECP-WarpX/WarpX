
#include "ImplicitSolverEM.H"
#include "WarpX.H"

void ImplicitSolverEM::Define ( WarpX* const  a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "ImplicitSolverEM object is already defined!");
    
    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;
 
    // Define E vectors
    m_E.Define( m_WarpX->Efield_fp );
    m_Eold.Define( m_WarpX->Efield_fp );
    SetDotMask(m_E);
    
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) { 
        // Define Bold vector
        const int lev = 0;
        m_Bold.resize(1); // size is number of levels
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& Bfp = *(m_WarpX->Bfield_fp[lev][n]);
            m_Bold[lev][n] = std::make_unique<amrex::MultiFab>( Bfp.boxArray(), 
                                                                Bfp.DistributionMap(),
                                                                Bfp.nComp(), 
                                                                Bfp.nGrowVect() );
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
        m_nlsolver = std::make_unique<PicardSolver<WarpXSolverVec,ImplicitSolverEM>>();
        m_nlsolver->Define(m_E, this);
    }
    else if (m_nlsolver_type == NonlinearSolverType::Newton) {
        m_nlsolver = std::make_unique<NewtonSolver<WarpXSolverVec,ImplicitSolverEM>>();
        m_nlsolver->Define(m_E, this);
    }

    m_is_defined = true;
}

void ImplicitSolverEM::PrintParameters () const
{
    if (!m_verbose) { return; }
    amrex::Print() << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "-------------- IMPLICIT EM SOLVER PARAMETERS --------------" << std::endl;
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

void ImplicitSolverEM::OneStep ( const amrex::Real  a_old_time,
                                 const amrex::Real  a_dt,
                                 const int          a_step )
{  
    using namespace amrex::literals;
    amrex::ignore_unused(a_step);

    // We have E^{n}.
    // Particles have p^{n} and x^{n}.
    // With theta implicit, B^{n}
    // With semi implicit, B^{n-1/2}

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart ( );

    // Save the fields at the start of the step
    m_Eold.Copy(m_WarpX->Efield_fp);
    m_E = m_Eold; // initial guess for E

    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        const int lev = 0;
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
    m_WarpX->UpdateElectricField( m_E, false ); // JRA not sure about false here. is what DG had.
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        m_WarpX->FinishMagneticField( m_Bold, m_theta );
    }
  
}

void ImplicitSolverEM::PreRHSOp ( const WarpXSolverVec&  a_E,
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

void ImplicitSolverEM::ComputeRHS ( WarpXSolverVec&  a_Erhs,
                              const WarpXSolverVec&  a_E,
                              const amrex::Real      a_time,
                              const amrex::Real      a_dt )
{  
    amrex::ignore_unused(a_E, a_time);
    m_WarpX->ComputeRHSE(m_theta*a_dt, a_Erhs);
}

void ImplicitSolverEM::UpdateWarpXState ( const WarpXSolverVec&  a_E,
                                          const amrex::Real      a_time,
                                          const amrex::Real      a_dt )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_time);

    // Update Efield_fp owned by WarpX
    m_WarpX->UpdateElectricField( a_E, true );

    // Update Bfield owned by WarpX
    if (m_WarpX->evolve_scheme == EvolveScheme::ThetaImplicit) {
        m_WarpX->UpdateMagneticField( m_Bold, m_theta*a_dt );
    }

    // The B field update needs. Talk to DG about this. Only needed when B updates?
    if (m_WarpX->num_mirrors>0){
        m_WarpX->applyMirrors(a_time);
        // E : guard cells are NOT up-to-date from the mirrors
        // B : guard cells are NOT up-to-date from the mirrors
    }

}

void ImplicitSolverEM::SetDotMask ( const WarpXSolverVec&  a_E ) {
    
    if (m_dot_mask_defined) { return; }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_E.IsDefined(),
        "ImplicitSolverEM::SetDotMask(a_E) called with undefined a_E");

    const amrex::Vector<amrex::Geometry>& Geom = m_WarpX->Geom();
    const int num_levels = a_E.getVec().size();
    m_dotMask.resize(num_levels);
    const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& solver_vec = a_E.getVec();
    for ( int n = 0; n < 3; n++) {
        const amrex::BoxArray& grids = solver_vec[0][n]->boxArray();
        amrex::MultiFab tmp( grids, solver_vec[0][n]->DistributionMap(),
                             1, 0, amrex::MFInfo().SetAlloc(false) );
        const amrex::Periodicity& period = Geom[0].periodicity();
        m_dotMask[0][n] = tmp.OwnerMask(period);
    }
    m_dot_mask_defined = true;

}

amrex::Real ImplicitSolverEM::dotProduct ( const WarpXSolverVec& a_X, 
                                           const WarpXSolverVec& a_Y ) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_dot_mask_defined,
        "ImplicitSolverEM::dotProduct called with m_dotMask not yet defined");
    amrex::Real result = 0.0;
    const int lev = 0;
    const bool local = true;
    for (int n = 0; n < 3; ++n) {
        auto rtmp = amrex::MultiFab::Dot( *m_dotMask[lev][n],
                                          *a_X.getVec()[lev][n], 0,
                                          *a_Y.getVec()[lev][n], 0, 1, 0, local);
        result += rtmp;
    }
    amrex::ParallelAllReduce::Sum(result, amrex::ParallelContext::CommunicatorSub());
    return result;
}

amrex::Real ImplicitSolverEM::norm ( const WarpXSolverVec& a_X ) const 
{
   amrex::Real result = dotProduct(a_X,a_X);
   return std::sqrt(result);
}
