#ifndef _IMPLICIT_SOLVER_EM_H_
#define _IMPLICIT_SOLVER_EM_H_

#include "FieldSolver/ImplicitSolvers/Utils/WarpXFieldVec.H"
#include "NonlinearSolvers/NonlinearSolverLibrary.H"

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_REAL.H>

class WarpX;
class ImplicitSolverEM
{
    public:

    ImplicitSolverEM() 
    {
        m_nlsolver_type = NonlinearSolverType::Picard;
        //m_nlsolver = nullptr;
        //m_func = nullptr;

        m_theta = 0.5;
        m_is_defined = false;
        m_semi_implicit = false;
        m_verbose = false;

        m_require_convergence = true;
        m_rtol = 1.0e-8;
        m_max_iter = 11;
    }

    ~ImplicitSolverEM() = default;
    
    //void Define();
    void Define( WarpX* const  a_WarpX );
    
    bool IsDefined() const { return m_is_defined; }
    
    amrex::Real theta() const { return m_theta; }
    
    void PrintParams() const;

    void Initialize();
    
    void OneStep( const amrex::Real  a_time,
                  const amrex::Real  a_dt,
                  const int          a_step );

    void PreRHSOp( const WarpXFieldVec&  a_Efield_vec,
                   const amrex::Real     a_time,
                   const amrex::Real     a_dt,
                   const int             a_nl_iter );

    void ComputeRHS( WarpXFieldVec&  a_Erhs_vec,
               const WarpXFieldVec&  a_Efield_vec,
               const amrex::Real     a_time,
               const amrex::Real     a_dt );

    void UpdateState( WarpXFieldVec&  a_Efield_vec,
                const amrex::Real     a_time );

    void PostUpdateState( const WarpXFieldVec&  a_field_vec,
                          const amrex::Real     a_time,
                          const amrex::Real     a_dt );

    private:

    bool m_is_defined;
    bool m_semi_implicit;     
    bool m_verbose;

    WarpX* m_WarpX;

    amrex::Real m_theta;
    
    NonlinearSolverType m_nlsolver_type;

    WarpXFieldVec m_E, m_Eold, m_Esave, m_Erhs;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_Bold;

    // picard solver parameters
    bool m_require_convergence;
    amrex::Real m_rtol;
    int m_max_iter;

};

#endif
