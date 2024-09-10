/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_LABFRAME_EXPLICITES_H_
#define WARPX_LABFRAME_EXPLICITES_H_

#include "ElectrostaticSolver.H"

class LabFrameExplicitES final : public ElectrostaticSolver
{
public:

    LabFrameExplicitES (int nlevs_max) : ElectrostaticSolver (nlevs_max) {ReadParameters();}

    void ComputeSpaceChargeField (
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_fp,
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_cp,
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& charge_buf,
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& phi_fp,
        [[maybe_unused]] MultiParticleContainer& mpc,
        [[maybe_unused]] MultiFluidContainer* mfl,
        [[maybe_unused]] amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Efield_fp,
        [[maybe_unused]] amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Bfield_fp
    ) override;

    void computePhiTriDiagonal (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi
    );

};

#endif
