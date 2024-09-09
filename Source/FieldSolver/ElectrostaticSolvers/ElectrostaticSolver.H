/* Copyright 2024 Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef Electrostatic_Solver_H_
#define Electrostatic_Solver_H_

#inlcude "ExplicitES.H"
#include "WarpX.H"

class ElectrostaticSolver
{
public:

    ElectrostaticSolver() = default;

    ComputeRHS() = default;

    ExplicitES m_explicit_es;

    ComputeSpaceChargeField (
        std::array<amrex::Vector<std::unique_ptr<amrex::MultiFab> >, 3>& E_field
    );

}

#endif
