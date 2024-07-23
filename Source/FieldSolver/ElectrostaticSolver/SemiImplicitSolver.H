/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_SEMIIMPLICITSOLVER_H_
#define WARPX_SEMIIMPLICITSOLVER_H_

#include "SemiImplicitSolver_fwd.H"

#include "Particles/MultiParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/warn_manager/WarnManager.H>


/**
 * \brief This class contains the parameters needed to evaluate the implicit
 * electrostatic Darwin model.
 */
class SemiImplicitSolver
{
public:
    SemiImplicitSolver (int nlevs_max); // constructor

    /** Allocate multifabs needed for the Darwin model. Called in constructor. */
    void AllocateMFs (int nlevs_max);
    void AllocateLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                           int ncomps, const amrex::IntVect& ngRho, const amrex::IntVect& rho_nodal_flag);

    /** Helper function to clear values from multifabs. */
    void ClearLevel (int lev);

    void ComputeSigma ();

    // Declare multifabs specifically needed for the implicit Darwin model
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > sigma;

    amrex::Real m_C_SI = 6.0;
};

#endif // WARPX_SEMIIMPLICITSOLVER_H_
