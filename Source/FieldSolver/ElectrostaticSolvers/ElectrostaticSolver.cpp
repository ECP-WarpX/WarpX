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
        m_explicit_es = std::make_unique<ExplicitES>(nlevs_max);
    }
    // else if (WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic)
    // {

    // }

    // Create an instance of the boundary handler to properly set boundary
    // conditions
    m_poisson_boundary_handler = std::make_unique<PoissonBoundaryHandler>();
}

ElectrostaticSolver::ComputeSpaceChargeField (
    std::array<amrex::Vector<std::unique_ptr<amrex::MultiFab> >, 3>& E_field,
)
{
    if (WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame)
    {
        m_explicit_es->ComputeSpaceChargeField(E_field);
    }
}
