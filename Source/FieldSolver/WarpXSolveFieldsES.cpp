/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Remi Lehe, Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "FieldSolver/ElectrostaticSolvers/ElectrostaticSolver.H"

void
WarpX::ComputeSpaceChargeField (bool const reset_fields)
{
    WARPX_PROFILE("WarpX::ComputeSpaceChargeField");
    if (reset_fields) {
        // Reset all E and B fields to 0, before calculating space-charge fields
        WARPX_PROFILE("WarpX::ComputeSpaceChargeField::reset_fields");
        for (int lev = 0; lev <= max_level; lev++) {
            for (int comp=0; comp<3; comp++) {
                Efield_fp[lev][comp]->setVal(0);
                Bfield_fp[lev][comp]->setVal(0);
            }
        }
    }

    m_electrostatic_solver->ComputeSpaceChargeField(
        rho_fp, rho_cp, charge_buf, phi_fp, *mypc, myfl.get(), Efield_fp, Bfield_fp
    );
}
