/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PICARD_SOLVER_H_
#define PICARD_SOLVER_H_

#include "NonlinearSolver.H"

#include <AMReX_ParmParse.H>
#include "Utils/TextMsg.H"

#include <vector>

/**
 * \brief Picard fixed-point iteration method to solve nonlinear
 *  equation of form: U = b + R(U)
 */

template<class Vec, class Ops>
class PicardSolver : public NonlinearSolver<Vec,Ops>
{
public:

    PicardSolver<Vec,Ops>() = default;

    ~PicardSolver<Vec,Ops>() override = default;

    // Prohibit Move and Copy operations
    PicardSolver(const PicardSolver&) = delete;
    PicardSolver& operator=(const PicardSolver&) = delete;
    PicardSolver(PicardSolver&&) noexcept = delete;
    PicardSolver& operator=(PicardSolver&&) noexcept = delete;

    void Define ( const Vec&  a_U,
                        Ops*  a_ops ) override;

    void Solve ( Vec&   a_U,
           const Vec&   a_b,
           amrex::Real  a_time,
           amrex::Real  a_dt ) const override;

    [[nodiscard]] bool IsDefined () const override { return m_is_defined; }

    void GetSolverParams ( amrex::Real&  a_rtol,
                           amrex::Real&  a_atol,
                           int&          a_maxits ) override
    {
        a_rtol = m_rtol;
        a_atol = m_atol;
        a_maxits = m_maxits;
    }

    void PrintParams () const override
    {
        amrex::Print() << "Picard max iterations:      " << m_maxits << std::endl;
        amrex::Print() << "Picard relative tolerance:  " << m_rtol << std::endl;
        amrex::Print() << "Picard absolute tolerance:  " << m_atol << std::endl;
        amrex::Print() << "Picard require convergence: " << (m_require_convergence?"true":"false") << std::endl;
    }

    void Verbose ( bool  a_verbose ) override { m_verbose = a_verbose; }

protected:

    bool m_is_defined = false;
    bool m_verbose = true;
    bool m_require_convergence = true;

    mutable Vec m_Usave, m_R;
    Ops* m_ops = nullptr;

/**
 * \brief relative tolerance for the Picard nonlinear solver
 */
    amrex::Real m_rtol = 1.0e-6;

/**
 * \brief absolute tolerance for the Picard nonlinear solver
 */
    amrex::Real m_atol = 0.;

/**
 * \brief maximum iterations for the Picard nonlinear solver
 */
    int m_maxits = 100;

private:

    void ParseParameters( );

};

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::Define ( const Vec&  a_U,
                                     Ops*        a_ops )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "Picard nonlinear solver object is already defined!");

    ParseParameters();

    m_Usave.Define(a_U);
    m_R.Define(a_U);

    m_ops = a_ops;

    m_is_defined = true;

}

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::ParseParameters ()
{
    const amrex::ParmParse pp_picard("picard");
    pp_picard.query("verbose",             m_verbose);
    pp_picard.query("absolute_tolerance",  m_atol);
    pp_picard.query("relative_tolerance",  m_rtol);
    pp_picard.query("max_iterations",      m_maxits);
    pp_picard.query("require_convergence", m_require_convergence);

}

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::Solve ( Vec&         a_U,
                              const Vec&         a_b,
                                    amrex::Real  a_time,
                                    amrex::Real  a_dt ) const
{
    BL_PROFILE("PicardSolver::Solve()");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_is_defined,
        "PicardSolver::Solve() called on undefined object");
    using namespace amrex::literals;

    //
    // Picard fixed-point iteration method to solve nonlinear
    // equation of form: U = b + R(U)
    //

    amrex::Real norm_abs = 0.;
    amrex::Real norm0 = 1._rt;
    amrex::Real norm_rel = 0.;

    int iter;
    for (iter = 0; iter < m_maxits;) {

        // Save previous state for norm calculation
        m_Usave = a_U;

        // Update the solver state (a_U = a_b + m_R)
        m_ops->PreRHSOp( a_U, a_time, a_dt, iter, false );
        m_ops->ComputeRHS( m_R, a_U, a_time, a_dt );
        a_U = a_b;
        a_U += m_R;

        // Compute the step norm and update iter
        m_Usave -= a_U;
        norm_abs = m_Usave.norm2();
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

        if (norm_abs < m_atol) {
            amrex::Print() << "Picard: exiting at iter = " << std::setw(3) << iter
                           << ". Satisfied absolute tolerance " << m_atol << std::endl;
            break;
        }

        if (norm_rel < m_rtol) {
            amrex::Print() << "Picard: exiting at iter = " << std::setw(3) << iter
                           << ". Satisfied relative tolerance " << m_rtol << std::endl;
            break;
        }

        if (iter >= m_maxits) {
            amrex::Print() << "Picard: exiting at iter = " << std::setw(3) << iter
                           << ". Maximum iteration reached: iter = " << m_maxits << std::endl;
            break;
        }

    }

    if (m_rtol > 0. && iter == m_maxits) {
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

}

#endif
