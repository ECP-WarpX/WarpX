/* Copyright 2024 Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ELECTROSTATICSOLVER_H_
#define WARPX_ELECTROSTATICSOLVER_H_

#include "ElectrostaticBase.H"
#include "LabFrameExplicitES.H"
#include "RelativisticExplicitES.H"

class ElectrostaticSolver
{
public:

    ElectrostaticSolver (int nlevs_max); // constructor

    void ClearLevel ( [[maybe_unused]] int lev ) {}
    void AllocateLevelMFs (
        // int lev, const BoxArray& ba, const DistributionMapping& dm,
        // const int ncomps, const IntVect& ngJ, const IntVect& rho_nodal_flag
    ) {}


    void ComputeSpaceChargeField (
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_fp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_cp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& charge_buf,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& phi_fp,
        MultiParticleContainer& mpc,
        MultiFluidContainer* mfl,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Efield_fp,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Bfield_fp
    );

    std::unique_ptr<ElectrostaticBase> m_electrostatic_solver;
};

#endif // WARPX_ELECTROSTATICSOLVER_H_
