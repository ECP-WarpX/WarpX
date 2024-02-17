#ifndef _NEWTON_SOLVER_H_
#define _NEWTON_SOLVER_H_

#include "NonlinearSolver.H"

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
            m_is_defined = false;
        }

        virtual ~NewtonSolver<Vec,Ops>() { }

        virtual void Define( const Vec&  a_U, 
                             Ops* const  a_ops );
                              //NonlinearFunction<Vec,Ops>* const )

        virtual void Solve( Vec&         a_U,
                      const Vec&         a_R,
                      const amrex::Real  a_time,
                      const amrex::Real  a_dt ) const;

        inline bool IsDefined() const { return m_is_defined; }

        inline void GetParams( amrex::Real&  a_rtol,
                               amrex::Real&  a_atol,
                               int&          a_maxits )
        {
            a_rtol = m_rtol;
            a_atol = m_atol;
            a_maxits = m_maxits;
        }
    
        void PrintParams( ) const
        {
            amrex::Print()     << "Newton absolute tolerance:  " << m_atol << std::endl;
            amrex::Print()     << "Newton relative tolerance:  " << m_rtol << std::endl;
            amrex::Print()     << "Newton max itsations:      " << m_maxits << std::endl;
            amrex::Print()     << "Newton require convergence: " << (m_require_convergence?"true":"false") << std::endl;
        }

        virtual void Verbose( const bool  a_verbose ) { m_verbose = a_verbose; }

    protected:

        bool m_is_defined;
        bool m_verbose;
        bool m_require_convergence;

        mutable Vec m_dU, m_F;
        mutable Ops* m_ops;

        amrex::Real m_atol, m_rtol;
        int    m_maxits;

    private:

        void ParseParameters( );

};

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::Define( const Vec&  a_U,
                                    Ops* const  a_ops )
                                    //NonlinearFunction<Vec,Ops>* const )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "Newton nonlinear solver object is already defined!");

    /* default */
    m_rtol = 1e-6;
    m_atol = 1e-12;
    m_maxits = 11;
    m_verbose = true;
    m_require_convergence = true;
  
    ParseParameters();
  
    m_dU.Define(a_U);
    m_F.Define(a_U);
  
    m_ops = a_ops;

    m_is_defined = true;

}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::ParseParameters( )
{
    const amrex::ParmParse pp_newton("newton");
    pp_newton.query("verbose",             m_verbose);
    pp_newton.query("absolute_tolerance",  m_atol);
    pp_newton.query("relative_tolerance",  m_rtol);
    pp_newton.query("max_itsations",      m_maxits);
    pp_newton.query("require_convergence", m_require_convergence);

}

template <class Vec, class Ops>
void NewtonSolver<Vec,Ops>::Solve( Vec&         a_U,
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
    // Newton routine to solve nonlinear
    // equation of form: m_F(a_U) = (a_U - G(a_U)) - a_b = 0
    //
    
    amrex::Real norm_abs, norm0, norm_rel;

    int its;
    for (its=0; its < m_maxits;) {

        //m_func->evalFunction(m_F, a_U, a_b, its);

        // Compute norm of the residual
        norm_abs = m_F.norm();
        if (its == 0) {
            if (norm_abs > 0.) { norm0 = norm_abs; }
            else { norm0 = 1._rt; }
        }
        norm_rel = norm_abs/norm0;
        its++;      
 
        // Check for convergence critsia
        if (m_verbose || its == m_maxits) {
            amrex::Print() << "Newton: its = " << std::setw(3) << its <<  ", norm = " 
                           << std::scientific << std::setprecision(5) << norm_abs << " (abs.), " 
                           << std::scientific << std::setprecision(5) << norm_rel << " (rel.)" << "\n";
        }

        if (norm_abs < m_rtol) {
            amrex::Print() << "Newton: exiting at its = " << std::setw(3) << its 
                           << ". Satisified absolute tolerance " << m_atol << std::endl;
            break;
        }
    
        if (norm_rel < m_rtol) {
            amrex::Print() << "Newton: exiting at its = " << std::setw(3) << its 
                           << ". Satisified relative tolerance " << m_rtol << std::endl;
            break;
        }
        
        //if (norm_abs > 100._rt*norm0) {
        //    amrex::Print() << "Newton: exiting at its = " << std::setw(3) << its 
        //                   << ". SOLVER DIVERGED " << m_rtol << std::endl;
        //    WARPX_ABORT_WITH_MESSAGE(convergenceMsg.str());
        //    break;
        //}

        // Solve linear system for Newton step [Jac]*dU = F
        m_dU.zero();
        //m_linear_solver->solve(m_dU, m_F);
    
        // Update solution U = U - dU
        a_U -= m_dU;

        if (its >= m_maxits) {
            amrex::Print() << "Newton: exiting at its = " << std::setw(3) << its 
                           << ". Maximum itsation reached: its = " << m_maxits << std::endl;
            break;
        }

    }

    if (m_rtol > 0. && m_atol > 0. && its == m_maxits) {
       std::stringstream convergenceMsg;
       convergenceMsg << "Newton solver failed to converge after " << its <<
                         " itsations. Relative norm is " << norm_rel <<
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

#endif
