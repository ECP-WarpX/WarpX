/* Copyright 2024 Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef Electrostatic_Solver_H_
#define Electrostatic_Solver_H_

#include "LabFrameExplicitES.H"

class ElectrostaticSolver
{
public:

    ElectrostaticSolver (int nlevs_max); // constructor

    std::unique_ptr<LabFrameExplicitES> m_labframe_explicit_es;

    void ComputeSpaceChargeField (
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > rho_fp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > rho_cp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > charge_buf,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > phi_fp,
        MultiParticleContainer& mpc,
        MultiFluidContainer* mfl,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> > Efield_fp,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> > Bfield_fp
    );
};

#endif
