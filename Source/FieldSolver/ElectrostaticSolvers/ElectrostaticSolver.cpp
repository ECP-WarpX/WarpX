/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ElectrostaticSolver.H"

ElectrostaticSolver::ElectrostaticSolver (int nlevs_max)
{
    if (WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame)
    {
        std::unique_ptr<LabFrameExplicitES> labframe_explicit_es = std::make_unique<LabFrameExplicitES>(nlevs_max);
        m_electrostatic_solver = std::move(labframe_explicit_es);
    }
    else if (WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic)
    {
        std::unique_ptr<RelativisticExplicitES> relativistic_explicit_es = std::make_unique<RelativisticExplicitES>(nlevs_max);
        m_electrostatic_solver = std::move(relativistic_explicit_es);
    }

}

void ElectrostaticSolver::ComputeSpaceChargeField (
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_fp,
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_cp,
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& charge_buf,
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& phi_fp,
    MultiParticleContainer& mpc,
    MultiFluidContainer* mfl,
    amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Efield_fp,
    amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Bfield_fp
)
{
    m_electrostatic_solver->ComputeSpaceChargeField(
        rho_fp, rho_cp, charge_buf, phi_fp, mpc, mfl, Efield_fp, Bfield_fp
    );
}
