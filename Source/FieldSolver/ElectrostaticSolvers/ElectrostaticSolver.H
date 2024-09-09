/* Copyright 2024 Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef Electrostatic_Solver_H_
#define Electrostatic_Solver_H_

class ElectrostaticSolver
{
    public:

    ElectrostaticSolver() = default;

    ComputeRHS() = default;

    PoissonBoundaryHandler m_poisson_boundary_handler;

    ExplicitES m_explicit_es;
}

#endif

