/* Copyright 2020
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#ifndef WARPX_DIM_RZ
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#else
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#endif
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <memory>

/**
 * Determine the timestep of the simulation. */
void
WarpX::ComputeDt ()
{
    // Determine
    const amrex::Real* dx = geom[max_level].CellSize();
    amrex::Real deltat = 0.;

    if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        // Computation of dt for spectral algorithm
        // (determined by the minimum cell size in all directions)
#if defined(WARPX_DIM_1D_Z)
        deltat = cfl * dx[0] / PhysConst::c;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        deltat = cfl * std::min(dx[0], dx[1]) / PhysConst::c;
#else
        deltat = cfl * std::min(dx[0], std::min(dx[1], dx[2])) / PhysConst::c;
#endif
    } else {
        // Computation of dt for FDTD algorithm
#ifdef WARPX_DIM_RZ
        // - In RZ geometry
        if (maxwell_solver_id == MaxwellSolverAlgo::Yee) {
            deltat = cfl * CylindricalYeeAlgorithm::ComputeMaxDt(dx,  n_rz_azimuthal_modes);
#else
        // - In Cartesian geometry
        if (do_nodal) {
            deltat = cfl * CartesianNodalAlgorithm::ComputeMaxDt(dx);
        } else if (maxwell_solver_id == MaxwellSolverAlgo::Yee
                    || maxwell_solver_id == MaxwellSolverAlgo::ECT) {
            deltat = cfl * CartesianYeeAlgorithm::ComputeMaxDt(dx);
        } else if (maxwell_solver_id == MaxwellSolverAlgo::CKC) {
            deltat = cfl * CartesianCKCAlgorithm::ComputeMaxDt(dx);
#endif
        } else {
            amrex::Abort(Utils::TextMsg::Err(
                "ComputeDt: Unknown algorithm"));
        }
    }

    dt.resize(0);
    dt.resize(max_level+1,deltat);

    if (do_subcycling) {
        for (int lev = max_level-1; lev >= 0; --lev) {
            dt[lev] = dt[lev+1] * refRatio(lev)[0];
        }
    }

    if (do_electrostatic != ElectrostaticSolverAlgo::None) {
        for (int lev=0; lev<=max_level; lev++) {
            dt[lev] = const_dt;
        }
    }
}

void
WarpX::PrintDtDxDyDz ()
{
    for (int lev=0; lev <= max_level; lev++) {
        const amrex::Real* dx_lev = geom[lev].CellSize();
        amrex::Print() << "Level " << lev << ": dt = " << dt[lev]
#if defined(WARPX_DIM_1D_Z)
                       << " ; dz = " << dx_lev[0] << '\n';
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                       << " ; dx = " << dx_lev[0]
                       << " ; dz = " << dx_lev[1] << '\n';
#elif defined(WARPX_DIM_3D)
                       << " ; dx = " << dx_lev[0]
                       << " ; dy = " << dx_lev[1]
                       << " ; dz = " << dx_lev[2] << '\n';
#endif
    }
}
