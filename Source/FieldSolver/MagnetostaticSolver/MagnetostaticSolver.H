/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: S. Eric Clark (LLNL), Roelof Groenewald (TAE)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MAGNETOSTATICSOLVER_H_
#define WARPX_MAGNETOSTATICSOLVER_H_

#include "MagnetostaticSolver_fwd.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>
#include <AMReX_REAL.H>

using namespace amrex;


/**
 * \brief This class contains functionality involved in solving the self
 * magnetic field of a plasma due to internal currents, i.e., the magnetostatic
 * solver. The solver uses the plasma current to calculate the vector potential
 * (in the Coulomb gauge) and evaluates B = ∇ x A.
 */
class MagnetostaticSolver
{
public:
    MagnetostaticSolver (int nlevs_max); // constructor

    /** Read user-defined model parameters. Called in constructor. */
    void ReadParameters ();

    /** Allocate member multifabs. Called in constructor. */
    void AllocateMFs (int nlevs_max);
    void AllocateLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                           int ncomps, const amrex::IntVect& ngJ, const amrex::IntVect& ngEB,
                           const amrex::IntVect& Ex_nodal_flag, const amrex::IntVect& Ey_nodal_flag,
                           const amrex::IntVect& Ez_nodal_flag);

    /** Helper function to clear values from member multifabs. */
    void ClearLevel (int lev);

    void InitData ();

    void UpdateEfromInductance (amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>& Efield,
                                amrex::Array<std::unique_ptr<amrex::MultiFab>, 3> const& Afield,
                                Real dt, const int lev);

    /** Boundary Handler for the Vector Potential Poisson Solver
     * This only will handle homogeneous Dirichlet boundary conditions on
     * embedded boundaries, and homogeneous dirichlet/Neumann or periodic
     * boundary conditions on domain boundaries.
     */
    struct BoundaryHandler {
        amrex::Array<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>, 3> lobc, hibc;
        bool bcs_set = false;
        std::array<std::array<bool, AMREX_SPACEDIM * 2>, 3> dirichlet_flag;
        bool has_non_periodic = false;

        void defineVectorPotentialBCs (
            Vector<FieldBoundaryType> field_boundary_lo,
            Vector<FieldBoundaryType> field_boundary_hi,
            [[maybe_unused]] Real prob_lo
        ) {
            for (int adim = 0; adim < 3; adim++) {
#ifdef WARPX_DIM_RZ
                int dim_start = 0;
                if (prob_lo == 0){
                    lobc[adim][0] = LinOpBCType::Neumann;
                    dirichlet_flag[adim][0] = false;
                    dim_start = 1;

                    // handle the r_max boundary explicitly
                    if (field_boundary_hi[0] == FieldBoundaryType::PEC) {
                        if (adim == 0) {
                            hibc[adim][0] = LinOpBCType::Neumann;
                            dirichlet_flag[adim][1] = false;
                        } else{
                            hibc[adim][0] = LinOpBCType::Dirichlet;
                            dirichlet_flag[adim][1] = true;
                        }
                    }
                    else if (field_boundary_hi[0] == FieldBoundaryType::Neumann) {
                        hibc[adim][0] = LinOpBCType::Neumann;
                        dirichlet_flag[adim][1] = false;
                    }
                }
#else
                const int dim_start = 0;
#endif
                bool ndotA = false;
                for (int idim=dim_start; idim<AMREX_SPACEDIM; idim++){
                    ndotA = (adim == idim);

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    if (idim == 1) { ndotA = (adim == 2); }
#endif

                    if ( field_boundary_lo[idim] == FieldBoundaryType::Periodic
                        && field_boundary_hi[idim] == FieldBoundaryType::Periodic ) {
                        lobc[adim][idim] = LinOpBCType::Periodic;
                        hibc[adim][idim] = LinOpBCType::Periodic;
                        dirichlet_flag[adim][idim*2] = false;
                        dirichlet_flag[adim][idim*2+1] = false;
                    }
                    else {
                        has_non_periodic = true;
                        if ( field_boundary_lo[idim] == FieldBoundaryType::PEC ) {
                            if (ndotA) {
                                lobc[adim][idim] = LinOpBCType::Neumann;
                                dirichlet_flag[adim][idim*2] = false;
                            } else {
                                lobc[adim][idim] = LinOpBCType::Dirichlet;
                                dirichlet_flag[adim][idim*2] = true;
                            }
                        }

                        else if ( field_boundary_lo[idim] == FieldBoundaryType::Neumann ) {
                            lobc[adim][idim] = LinOpBCType::Neumann;
                            dirichlet_flag[adim][idim*2] = false;
                        }
                        else {
                            WARPX_ABORT_WITH_MESSAGE(
                                "Field boundary conditions have to be either periodic, PEC or neumann "
                                "when using the magnetostatic solver."
                            );
                        }

                        if ( field_boundary_hi[idim] == FieldBoundaryType::PEC ) {
                            if (ndotA) {
                                hibc[adim][idim] = LinOpBCType::Neumann;
                                dirichlet_flag[adim][idim*2+1] = false;
                            } else {
                                hibc[adim][idim] = LinOpBCType::Dirichlet;
                                dirichlet_flag[adim][idim*2+1] = true;
                            }
                        }
                        else if ( field_boundary_hi[idim] == FieldBoundaryType::Neumann ) {
                            hibc[adim][idim] = LinOpBCType::Neumann;
                            dirichlet_flag[adim][idim*2+1] = false;
                        }
                        else {
                            WARPX_ABORT_WITH_MESSAGE(
                                "Field boundary conditions have to be either periodic, PEC or neumann "
                                "when using the magnetostatic solver."
                            );
                        }
                    }
                }
            }
            bcs_set = true;
        }
    };

    BoundaryHandler m_boundary_handler = BoundaryHandler();

    // numerical parameters
    amrex::Real required_precision = 1.e-11_rt;
    amrex::Real absolute_tolerance = 0.0_rt;
    int max_iters = 200;
    int verbosity = 2;

    // Declare multifabs specifically needed for the magnetostatic solver
    Vector<std::array< std::unique_ptr<MultiFab>, 3 > > Afield_fp_nodal;
    Vector<std::array< std::unique_ptr<MultiFab>, 3 > > Afield_fp_old;
    Vector<std::array< std::unique_ptr<MultiFab>, 3 > > current_fp_temp;
};
#endif //WARPX_MAGNETOSTATICSOLVER_H_
