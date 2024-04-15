/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_IMPLICITDARWINSOLVER_H_
#define WARPX_IMPLICITDARWINSOLVER_H_

#include "ImplicitDarwinSolver_fwd.H"

#include "Particles/MultiParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>
#include <ablastr/fields/PoissonSolver.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>

/**
 * \brief This class contains the parameters needed to evaluate the implicit
 * electrostatic Darwin model.
 */
class ImplicitDarwinSolver
{
public:
    ImplicitDarwinSolver (int nlevs_max); // constructor

    /** Allocate multifabs needed for the Darwin model. Called in constructor. */
    void AllocateMFs (int nlevs_max);
    void AllocateLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                           int ncomps, const amrex::IntVect& ngRho, const amrex::IntVect& rho_nodal_flag);

    /** Helper function to clear values from multifabs. */
    void ClearLevel (int lev);

    /**
     * \brief
     * Function to update the E-field using the implicit Darwin model.
     */
    void AddSpaceChargeField (
        amrex::Vector<std::unique_ptr<amrex::MultiFab>>& rhofield,
        amrex::Vector<std::unique_ptr<amrex::MultiFab>>& phifield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
        amrex::Real required_precision, amrex::Real absolute_tolerance,
        int max_iters, int verbosity
    );

    template<typename T_BoundaryHandler, typename T_FArrayBoxFactory = void>
    void ComputePhi (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
        amrex::Real relative_tolerance, amrex::Real absolute_tolerance,
        int max_iters, int verbosity,
        amrex::Geometry const geom,
        T_BoundaryHandler const boundary_handler,
        std::optional<amrex::Vector<T_FArrayBoxFactory const *> > eb_farray_box_factory = std::nullopt
    ) const;

    void ComputeSigma ();

    // Declare multifabs specifically needed for the implicit Darwin model
    // amrex::Vector<std::unique_ptr<amrex::MultiFab>> sigma_fp;
    // amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM > > sigma;
    // amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > > sigma;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > sigma;

    amrex::Real m_C_SI = 6.0;
};

#endif // WARPX_IMPLICITDARWINSOLVER_H_
