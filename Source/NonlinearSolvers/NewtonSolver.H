#ifndef _NEWTON_SOLVER_H_
#define _NEWTON_SOLVER_H_

#include "NonlinearSolver.H"
#include "JacobianFunctionJFNK.H"

#include <AMReX_GMRES.H>
#include <AMReX_ParmParse.H>
#include "Utils/TextMsg.H"

#include <vector>

template<class Vec, class Ops>
class NewtonSolver : public NonlinearSolver<Vec,Ops>
{
public:

    NewtonSolver<Vec,Ops>() 
    {
        m_ops = nullptr;
        m_update_pc = true;
        m_update_pc_init = true;
        m_is_defined = false;
    }

    virtual ~NewtonSolver<Vec,Ops>() { }

    virtual void Define ( const Vec&  a_U, 
                          Ops* const  a_ops );

    virtual void Solve ( Vec&         a_U,
                   const Vec&         a_R,
                   const amrex::Real  a_time,
                   const amrex::Real  a_dt ) const;

    inline bool IsDefined () const { return m_is_defined; }

    inline void GetParams ( amrex::Real&  a_rtol,
                            amrex::Real&  a_atol,
                            int&          a_maxits )
    {
        a_rtol = m_rtol;
        a_atol = m_atol;
        a_maxits = m_maxits;
    }
    
    inline void CurTime ( const amrex::Real  a_time ) const
    {
        m_cur_time = a_time;
        m_linear_function->curTime( a_time );
    }

    inline void CurTimeStep ( const amrex::Real  a_dt ) const
    {
        m_dt = a_dt;
        m_linear_function->curTimeStep( a_dt );
    }
    
    void PrintParams () const
    {
        amrex::Print()     << "Newton verbose:             " << (m_verbose?"true":"false") << std::endl;
        amrex::Print()     << "Newton max iterations:      " << m_maxits << std::endl;
        amrex::Print()     << "Newton relative tolerance:  " << m_rtol << std::endl;
        amrex::Print()     << "Newton absolute tolerance:  " << m_atol << std::endl;
        amrex::Print()     << "Newton require convergence: " << (m_require_convergence?"true":"false") << std::endl;
        amrex::Print()     << "GMRES verbose:            " << m_gmres_verbose_int << std::endl;
        amrex::Print()     << "GMRES max iterations:     " << m_gmres_maxits << std::endl;
        amrex::Print()     << "GMRES relative tolerance: " << m_gmres_rtol << std::endl;
        amrex::Print()     << "GMRES absolute tolerance: " << m_gmres_atol << std::endl;
    }

    virtual void Verbose ( const bool  a_verbose ) { m_verbose = a_verbose; }

protected:

    bool m_is_defined;
    bool m_verbose;
    bool m_require_convergence;

    mutable Vec m_dU, m_F, m_R;
    Ops* m_ops;

    amrex::Real m_atol, m_rtol;
    int m_maxits;
    
    amrex::Real m_gmres_atol, m_gmres_rtol;
    int m_gmres_verbose_int, m_gmres_maxits;
    
    mutable amrex::Real m_cur_time, m_dt;
    mutable bool m_update_pc;
    mutable bool m_update_pc_init;

private:

    std::unique_ptr<JacobianFunctionJFNK<Vec,Ops>> m_linear_function;
    std::unique_ptr<amrex::GMRES<Vec,JacobianFunctionJFNK<Vec,Ops>>> m_linear_solver;

    void ParseParameters ();

    void EvalResidual ( Vec&         a_R, 
                  const Vec&         a_U,
                  const Vec&         a_b,
                  const amrex::Real  a_time,
                  const amrex::Real  a_dt,
                  const int          a_iter ) const;

};

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::Define ( const Vec&  a_U,
                                     Ops* const  a_ops )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "Newton nonlinear solver object is already defined!");

    /* default */
    m_rtol = 1.0e-6;
    m_atol = 0.0;
    m_maxits = 100;
    m_verbose = true;
    m_require_convergence = true;
    
    m_gmres_verbose_int = 0;
    m_gmres_rtol = 1.0e-4;
    m_gmres_atol = 0.0;
    m_gmres_maxits = 1000;
  
    ParseParameters();
  
    m_dU.Define(a_U);
    m_F.Define(a_U); // residual function F(U) = U - b - R(U) = 0
    m_R.Define(a_U);
  
    m_ops = a_ops;

    m_linear_function = std::make_unique<JacobianFunctionJFNK<Vec,Ops>>();
    m_linear_function->define(m_F, m_ops);
    
    m_linear_solver = std::make_unique<amrex::GMRES<Vec,JacobianFunctionJFNK<Vec,Ops>>>();
    m_linear_solver->define(*m_linear_function);
    m_linear_solver->setVerbose( m_gmres_verbose_int );

    //
    // JRA, below is used to test GMRES, and linOp functions defined in JacobianFunctionJFNK
    //
    /*
    m_linear_solver->solve(m_dU,m_F,m_gmres_rtol,m_gmres_atol,m_gmres_maxits);
    
    JacobianFunctionJFNK<Vec,Ops> linOp;
    linOp.define(m_F,m_ops);
    m_F.setVal(1.0);
    m_dU.setVal(1.0);
    //const amrex::Real dotProd = linOp.dotProduct(m_F,m_dU);
    const amrex::Real dotProd = m_ops->dotProduct(m_F,m_dU);
    amrex::Print() << " JRA: dotProduct = " << dotProd << std::endl;
    m_dU.setVal(1.0);
    linOp.increment(m_F,m_dU,2.0);
    const amrex::Real dotProd2 = linOp.dotProduct(m_F,m_dU);
    amrex::Print() << " JRA: dotProduct2 = " << dotProd2 << std::endl;

    int gmres_iter, gmres_stat;
    gmres_iter = m_linear_solver->getNumIters();
    gmres_stat = m_linear_solver->getStatus();
    amrex::Print() << "JRA: gmres_iter = " << gmres_iter << std::endl;
    amrex::Print() << "JRA: gmres_stat = " << gmres_stat << std::endl;
    */
    //
    //
    //

    m_is_defined = true;

}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::ParseParameters ()
{
    const amrex::ParmParse pp_newton("newton");
    pp_newton.query("verbose",             m_verbose);
    pp_newton.query("absolute_tolerance",  m_atol);
    pp_newton.query("relative_tolerance",  m_rtol);
    pp_newton.query("max_iterations",      m_maxits);
    pp_newton.query("require_convergence", m_require_convergence);
    
    const amrex::ParmParse pp_gmres("gmres");
    pp_gmres.query("verbose_int",         m_gmres_verbose_int);
    pp_gmres.query("absolute_tolerance",  m_gmres_atol);
    pp_gmres.query("relative_tolerance",  m_gmres_rtol);
    pp_gmres.query("max_iterations",      m_gmres_maxits);
}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::Solve ( Vec&         a_U,
                              const Vec&         a_b,
                              const amrex::Real  a_time,
                              const amrex::Real  a_dt ) const
{
    BL_PROFILE("NewtonSolver::Solve()");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_is_defined,
        "NewtonSolver::Solve() called on undefined object");
    using namespace amrex::literals;

    //
    // Newton routine to solve nonlinear equation of form:
    // F(U) = U - b - R(U) = 0
    //
    
    CurTime(a_time);
    CurTimeStep(a_dt);

    amrex::Real norm_abs, norm0, norm_rel;
    int its;
    for (its=0; its < m_maxits;) {

        EvalResidual(m_F, a_U, a_b, a_time, a_dt, its);
        
        // Compute norm of the residual
        norm_abs = m_ops->norm(m_F);
        if (its == 0) {
            if (norm_abs > 0.) { norm0 = norm_abs; }
            else { norm0 = 1._rt; }
        }
        norm_rel = norm_abs/norm0;
 
        // Check for convergence critsia
        if (m_verbose || its == m_maxits) {
            amrex::Print() << "Newton: iteration = " << std::setw(3) << its <<  ", norm = " 
                           << std::scientific << std::setprecision(5) << norm_abs << " (abs.), " 
                           << std::scientific << std::setprecision(5) << norm_rel << " (rel.)" << "\n";
        }

        if (norm_abs < m_rtol) {
            amrex::Print() << "Newton: exiting at iteration = " << std::setw(3) << its 
                           << ". Satisfied absolute tolerance " << m_atol << std::endl;
            break;
        }
    
        if (norm_rel < m_rtol) {
            amrex::Print() << "Newton: exiting at iteration = " << std::setw(3) << its 
                           << ". Satisfied relative tolerance " << m_rtol << std::endl;
            break;
        }
        
        if (norm_abs > 100._rt*norm0) {
            amrex::Print() << "Newton: exiting at iteration = " << std::setw(3) << its 
                           << ". SOLVER DIVERGED! relative tolerance = " << m_rtol << std::endl;
            std::stringstream convergenceMsg;
            convergenceMsg << "Newton: exiting at iteration " << std::setw(3) << its <<
                              ". SOLVER DIVERGED! absolute norm = " << norm_abs <<
                              " has increased by 100X from that after first iteration.";
            WARPX_ABORT_WITH_MESSAGE(convergenceMsg.str());
        }

        // Solve linear system for Newton step [Jac]*dU = F
        m_dU.zero();
        m_linear_solver->solve(m_dU,m_F,m_gmres_rtol,m_gmres_atol,m_gmres_maxits);
    
        // Update solution
        a_U -= m_dU;

        its++;
        if (its >= m_maxits) {
            amrex::Print() << "Newton: exiting at its = " << std::setw(3) << its 
                           << ". Maximum iteration reached: its = " << m_maxits << std::endl;
            break;
        }

    }

    if (m_rtol > 0. && its == m_maxits) {
       std::stringstream convergenceMsg;
       convergenceMsg << "Newton solver failed to converge after " << its <<
                         " iterations. Relative norm is " << norm_rel <<
                         " and the relative tolerance is " << m_rtol <<
                         ". Absolute norm is " << norm_abs << 
                         " and the absolute tolerance is " << m_atol;
       if (m_verbose) { amrex::Print() << convergenceMsg.str() << std::endl; }
       if (m_require_convergence) {
           WARPX_ABORT_WITH_MESSAGE(convergenceMsg.str());
       } else {
           ablastr::warn_manager::WMRecordWarning("NewtonSolver", convergenceMsg.str());
       }
    }
 
    return;

}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::EvalResidual ( Vec&         a_F,
                                     const Vec&         a_U,
                                     const Vec&         a_b,
                                     const amrex::Real  a_time,
                                     const amrex::Real  a_dt,
                                     const int          a_iter ) const
{

    m_ops->PreRHSOp( a_U, a_time, a_dt, a_iter, false );
    m_ops->ComputeRHS( m_R, a_U, a_time, a_dt );

    // set base Y and R(Y) for Jacobian //
    m_linear_function->setBaseSolution(a_U);
    m_linear_function->setBaseRHS(m_R);
  
    // update preconditioner //
    if (m_update_pc || m_update_pc_init) {
        m_linear_function->updatePreCondMat(a_U);
    }
    m_update_pc_init = false;

    // Compute ([Y - R(Y)] - b) //
    a_F = a_U - m_R - a_b;

}

#endif
