#ifndef _PICARD_SOLVER_H_
#define _PICARD_SOLVER_H_

#include "NonlinearSolver.H"

#include <AMReX_ParmParse.H>
#include "Utils/TextMsg.H"

#include <vector>

template<class Vec, class Ops>
class PicardSolver : public NonlinearSolver<Vec,Ops>
{
    public:

        PicardSolver<Vec,Ops>() 
        {
            m_ops = nullptr;
            m_is_defined = false;
        }

        virtual ~PicardSolver<Vec,Ops>() { }

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
            amrex::Print()     << "Picard absolute tolerance:  " << m_atol << std::endl;
            amrex::Print()     << "Picard relative tolerance:  " << m_rtol << std::endl;
            amrex::Print()     << "Picard max iterations:      " << m_maxits << std::endl;
            amrex::Print()     << "Picard require convergence: " << (m_require_convergence?"true":"false") << std::endl;
        }

        virtual void Verbose( const bool  a_verbose ) { m_verbose = a_verbose; }

    protected:

        bool m_is_defined;
        bool m_verbose;
        bool m_require_convergence;

        mutable Vec m_Usave, m_F;
        mutable Ops* m_ops;

        amrex::Real m_atol, m_rtol;
        int    m_maxits;

    private:

        void ParseParameters( );

};

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::Define( const Vec&  a_U,
                                    Ops* const  a_ops )
                                    //NonlinearFunction<Vec,Ops>* const )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "Picard nonlinear solver object is already defined!");

    /* default */
    m_rtol = 1e-6;
    m_atol = 1e-12;
    m_maxits = 11;
    m_verbose = true;
    m_require_convergence = true;
  
    ParseParameters();
  
    m_Usave.Define(a_U);
    m_F.Define(a_U);
  
    m_ops = a_ops;
    
    m_is_defined = true;

}

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::ParseParameters( )
{
    const amrex::ParmParse pp_picard("picard");
    pp_picard.query("verbose",             m_verbose);
    pp_picard.query("absolute_tolerance",  m_atol);
    pp_picard.query("relative_tolerance",  m_rtol);
    pp_picard.query("max_iterations",      m_maxits);
    pp_picard.query("require_convergence", m_require_convergence);

}

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::Solve( Vec&         a_U,
                             const Vec&         a_R,
                             const amrex::Real  a_time,
                             const amrex::Real  a_dt ) const
{
    BL_PROFILE("PicardSolver::Solve()");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_is_defined,
        "PicardSolver::Solve() called on undefined object");
    using namespace amrex::literals;

    //
    // Picard routine to solve nonlinear
    // equation of form: a_U = a_R + m_F(a_U)
    //

    amrex::Real norm_abs, norm0, norm_rel;

    int iter(0);
    while(1) {

        // Save previous state for norm calculation
        m_Usave = a_U;

        // Update the solver state
        m_ops->PreRHSOp( a_U, a_time, a_dt, iter, false );
        m_ops->ComputeRHS( m_F, a_U, a_time, a_dt );
        a_U = a_R + m_F;
        
        // Do post update operations
        m_ops->PostUpdateState( a_U, a_time, a_dt );

        // Compute the step norm and update iter
        m_Usave -= a_U;
        //norm_abs = m_Usave.norm();
        norm_abs = m_ops->norm(m_Usave);
        if (iter == 0) {
            if (norm_abs > 0.) { norm0 = norm_abs; }
            else { norm0 = 1._rt; }
        }
        norm_rel = norm_abs/norm0;
        iter++;
        
        // Check for convergence criteria
        if (m_verbose || iter == m_maxits) {
            amrex::Print() << "Picard: iter = " << std::setw(3) << iter <<  ", norm = " 
                           << std::scientific << std::setprecision(5) << norm_abs << " (abs.), " 
                           << std::scientific << std::setprecision(5) << norm_rel << " (rel.)" << "\n";
        }

        if (norm_abs < m_rtol) {
            amrex::Print() << "Picard: exiting at iter = " << std::setw(3) << iter 
                           << ". Satisified absolute tolerance " << m_atol << std::endl;
            break;
        }
    
        if (norm_rel < m_rtol) {
            amrex::Print() << "Picard: exiting at iter = " << std::setw(3) << iter 
                           << ". Satisified relative tolerance " << m_rtol << std::endl;
            break;
        }

        if (iter >= m_maxits) {
            amrex::Print() << "Picard: exiting at iter = " << std::setw(3) << iter 
                           << ". Maximum iteration reached: iter = " << m_maxits << std::endl;
            break;
        }

    }

    if (m_rtol > 0. && m_atol > 0. && iter == m_maxits) {
       std::stringstream convergenceMsg;
       convergenceMsg << "Picard solver failed to converge after " << iter <<
                         " iterations. Relative norm is " << norm_rel <<
                         " and the relative tolerance is " << m_rtol <<
                         ". Absolute norm is " << norm_abs << 
                         " and the absolute tolerance is " << m_atol;
       if (m_verbose) { amrex::Print() << convergenceMsg.str() << std::endl; }
       if (m_require_convergence) {
           WARPX_ABORT_WITH_MESSAGE(convergenceMsg.str());
       } else {
           ablastr::warn_manager::WMRecordWarning("PicardSolver", convergenceMsg.str());
       }
    }

    return;

}

#endif
