/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SemiImplicitEM.H"
#include "WarpX.H"

/*
 *  Electromagnetic semi-implicit time solver class. The electric field and the
 *  particles are implicitly coupled in this algorithm, but the magnetic field
 *  is advanced in the standard explicit leap-frog manner (hence semi-implicit).
 *
 *  The time stencil is
 *  Eg^{n+1} = Eg^{n} + c^2*dt*( curlBg^{n+1/2} - mu0*Jg^{n+1/2} )
 *  Bg^{n+3/2} = Bg^{n+1/2} - dt*curlEg^{n+1}
 *  xp^{n+1} = xp^n + dt*up^{n+1/2}/(gammap^n + gammap^{n+1})
 *  up^{n+1} = up^n + dt*qp/mp*(Ep^{n+1/2} + up^{n+1/2}/gammap^{n+1/2} x Bp^{n+1/2})
 *  where f^{n+1/2} = (f^{n} + f^{n+1})/2.0, for all but Bg, which lives at half steps
 *
 *  This algorithm is approximately energy conserving. The violation in energy conservation
 *  is typically negligible. The advantage of this method over the exactly energy-conserving
 *  theta-implicit EM method is that light wave dispersion is captured much better. However,
 *  the CFL condition for light waves does have to be satisifed for numerical stability.
 *
 *  See G. Chen, L. Chacon, L. Yin, B.J. Albright, D.J. Stark, R.F. Bird,
 *  "A semi-implicit energy- and charge-conserving particle-in-cell algorithm for the
 *  relativistic Vlasov-Maxwell equations.", JCP 407 (2020).
 */

void SemiImplicitEM::Define ( WarpX*  a_WarpX )
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
    m_WarpX->ApplyMagneticFieldBCs();

    // Solve nonlinear system for E at t_{n+1/2}
    // Particles will be advanced to t_{n+1/2}
    m_nlsolver->Solve( m_E, m_Eold, a_time, a_dt );

    // update WarpX owned Efield_fp and Bfield_fp to t_{n+1/2}
    UpdateWarpXState( m_E, a_time, a_dt );

    // Update field boundary probes prior to updating fields to t_{n+1}
    //m_fields->updateBoundaryProbes( a_dt );

    // Advance particles to step n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance fields to step n+1
    m_WarpX->FinishImplicitField(m_E.getVec(), m_Eold.getVec(), 0.5_rt);
    m_WarpX->SetElectricFieldAndApplyBCs( m_E );

}

void SemiImplicitEM::PreRHSOp ( const WarpXSolverVec&  a_E,
                                amrex::Real            a_time,
                                amrex::Real            a_dt,
                                int                    a_nl_iter,
                                bool                   a_from_jacobian )
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
                                  amrex::Real      a_time,
                                  amrex::Real      a_dt )
{
    amrex::ignore_unused(a_E, a_time);
    using namespace amrex::literals;
    m_WarpX->ComputeRHSE(0.5_rt*a_dt, a_Erhs);
}

void SemiImplicitEM::UpdateWarpXState ( const WarpXSolverVec&  a_E,
                                        amrex::Real            a_time,
                                        amrex::Real            a_dt )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_time,a_dt);

    // Update Efield_fp owned by WarpX
    m_WarpX->SetElectricFieldAndApplyBCs( a_E );

    // The B field update needs. Talk to DG about this. Only needed when B updates?
    if (WarpX::num_mirrors>0){
        m_WarpX->applyMirrors(a_time);
        // E : guard cells are NOT up-to-date from the mirrors
        // B : guard cells are NOT up-to-date from the mirrors
    }

}
