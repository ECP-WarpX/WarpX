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

#include "ElectrostaticBase.H"

class LabFrameExplicitES : public ElectrostaticBase
{
public:

    void ReadParameters ();
    void ComputeSpaceChargeField (
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > rho_fp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > rho_cp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > charge_buf,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> > phi_fp,
        MultiParticleContainer& mpc,
        MultiFluidContainer* mfl,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> > Efield_fp
    );
    void computePhiTriDiagonal (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi
    );
    // Multifab for electrostatic potential
    // amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_phi_fp;
    // amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_rho_fp;
};

#endif
