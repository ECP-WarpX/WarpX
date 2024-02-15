/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

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

        NonlinearSolver<Vec,Ops>() { }
        virtual ~NonlinearSolver<Vec,Ops>() { }

        virtual void define(  const Vec&,
                              Ops* const,
                              //NonlinearFunction<Vec,Ops>* const,
                              const amrex::Real ) = 0;

        virtual void solve( Vec&,
                            const Vec&,
                            const amrex::Real,
                            const amrex::Real ) const = 0;

        virtual bool isDefined() const = 0;
        virtual void printParams() const = 0;
        virtual void getParams(double&, double&, int&) = 0;

        virtual void verbose( const bool ) { }
        virtual void setOutputIndent( const std::string& ) { }

        virtual void getConvergenceParams(  int& a_l_exit_status, int& a_l_last_iter, int& a_l_total_iter,
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

#endif // WARPX_NONLINEAR_SOLVER_H_
