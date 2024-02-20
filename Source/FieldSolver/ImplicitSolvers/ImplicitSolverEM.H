#ifndef _IMPLICIT_SOLVER_EM_H_
#define _IMPLICIT_SOLVER_EM_H_

#include "FieldSolver/ImplicitSolvers/Utils/WarpXSolverVec.H"
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
        //m_func = nullptr;

        m_theta = 0.5;
        m_is_defined = false;
        m_semi_implicit = false;
        m_verbose = false;

    }

    ~ImplicitSolverEM() = default;
    
    //void Define();
    void Define( WarpX* const  a_WarpX );
    
    bool IsDefined() const { return m_is_defined; }
    
    amrex::Real theta() const { return m_theta; }
    
    void PrintParameters() const;

    void Initialize();
    
    void OneStep( const amrex::Real  a_time,
                  const amrex::Real  a_dt,
                  const int          a_step );

    void PreRHSOp( const WarpXSolverVec&  a_E,
                   const amrex::Real      a_time,
                   const amrex::Real      a_dt,
                   const int              a_nl_iter );

    void ComputeRHS( WarpXSolverVec&  a_Erhs,
               const WarpXSolverVec&  a_E,
               const amrex::Real      a_time,
               const amrex::Real      a_dt );

    void PostUpdateState( const WarpXSolverVec&  a_E,
                          const amrex::Real      a_time,
                          const amrex::Real      a_dt );

    const amrex::Vector<amrex::Geometry>& getGeomVec() const { return m_geom_vec; }

    private:

    bool m_is_defined;
    bool m_semi_implicit;     
    bool m_verbose;

    WarpX* m_WarpX;

    amrex::Real m_theta;
    amrex::Vector<amrex::Geometry> m_geom_vec;

    NonlinearSolverType m_nlsolver_type;
    std::unique_ptr<NonlinearSolver<WarpXSolverVec,ImplicitSolverEM>> m_nlsolver;

    WarpXSolverVec m_E, m_Eold;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_Bold;

};

#endif
