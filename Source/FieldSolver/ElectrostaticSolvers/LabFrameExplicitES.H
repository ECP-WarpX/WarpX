/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_LABFRAMEEXPLICITES_H_
#define WARPX_LABFRAMEEXPLICITES_H_

#include "ElectrostaticSolver.H"

class LabFrameExplicitES final : public ElectrostaticSolver
{
public:

    LabFrameExplicitES (int nlevs_max) : ElectrostaticSolver (nlevs_max) {
        ReadParameters();
    }

    void InitData () override;

    void ComputeSpaceChargeField (
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_fp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_cp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& charge_buf,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& phi_fp,
        MultiParticleContainer& mpc,
        MultiFluidContainer* mfl,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Efield_fp,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Bfield_fp
    ) override;

    void computePhiTriDiagonal (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi
    );

};

#endif  // WARPX_LABFRAMEEXPLICITES_H_
