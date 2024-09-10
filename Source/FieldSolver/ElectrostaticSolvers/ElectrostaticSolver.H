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

    void ClearLevel ( int lev );
    void AllocateLevelMFs (
        int lev, const BoxArray& ba, const DistributionMapping& dm,
        const int ncomps, const IntVect& ngJ, const IntVect& rho_nodal_flag
    );


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

    std::unique_ptr<LabFrameExplicitES> m_labframe_explicit_es;
    std::unique_ptr<RelativisticExplicitES> m_relativistic_explicit_es;
};

#endif // WARPX_ELECTROSTATICSOLVER_H_
