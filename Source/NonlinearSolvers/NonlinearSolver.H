#ifndef WARPX_NONLINEAR_SOLVER_H_
#define WARPX_NONLINEAR_SOLVER_H_

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <string>
#include <array>
#include <memory>

/**
 * \brief Top-level class for the nonlinear solver
 */
template<class Vec, class Ops>
class NonlinearSolver
{
    public:

        NonlinearSolver<Vec,Ops>() = default;
        virtual ~NonlinearSolver<Vec,Ops>() = default;

        virtual void Define ( const Vec&,
                              Ops* ) = 0;

        virtual void Solve ( Vec&,
                       const Vec&,
                             amrex::Real,
                             amrex::Real ) const = 0;

        virtual bool IsDefined () const = 0;
        virtual void PrintParams () const = 0;
        virtual void GetParams (amrex::Real&, amrex::Real&, int&) = 0;

        virtual void Verbose ( bool ) { }

        virtual void GetConvergenceParams (  int& a_l_exit_status, int& a_l_last_iter, int& a_l_total_iter,
                                             int& a_nl_exit_status, int& a_nl_iter, int& a_nl_total_iter,
                                          double& a_nl_abs_res, double& a_nl_rel_res ) const
        {
            a_l_exit_status = 0;
            a_l_last_iter = 0;
            a_l_total_iter = 0;
            a_nl_exit_status = 0;
            a_nl_iter = 0;
            a_nl_total_iter = 0;
            a_nl_abs_res = 0.0;
            a_nl_rel_res = 0.0;
        }

    private:

};

#endif
