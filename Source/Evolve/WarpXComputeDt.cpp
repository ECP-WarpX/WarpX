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
    // Handle cases where the timestep is not limited by the speed of light
    if (electromagnetic_solver_id == ElectromagneticSolverAlgo::None) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_const_dt.has_value(), "warpx.const_dt must be specified with the electrostatic solver.");
        for (int lev=0; lev<=max_level; lev++) {
            dt[lev] = m_const_dt.value();
        }
        return;
    }

    // Determine the appropriate timestep as limited by the speed of light
    const amrex::Real* dx = geom[max_level].CellSize();
    amrex::Real deltat = 0.;

    if (m_const_dt.has_value()) {
        deltat = m_const_dt.value();
    } else if (electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD) {
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
        if (electromagnetic_solver_id == ElectromagneticSolverAlgo::Yee) {
            deltat = cfl * CylindricalYeeAlgorithm::ComputeMaxDt(dx,  n_rz_azimuthal_modes);
#else
        // - In Cartesian geometry
        if (do_nodal) {
            deltat = cfl * CartesianNodalAlgorithm::ComputeMaxDt(dx);
        } else if (electromagnetic_solver_id == ElectromagneticSolverAlgo::Yee
                    || electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT) {
            deltat = cfl * CartesianYeeAlgorithm::ComputeMaxDt(dx);
        } else if (electromagnetic_solver_id == ElectromagneticSolverAlgo::CKC) {
            deltat = cfl * CartesianCKCAlgorithm::ComputeMaxDt(dx);
#endif
        } else {
            amrex::Abort(Utils::TextMsg::Err("ComputeDt: Unknown algorithm"));
        }
    }

    dt.resize(0);
    dt.resize(max_level+1,deltat);

    if (do_subcycling) {
        for (int lev = max_level-1; lev >= 0; --lev) {
            dt[lev] = dt[lev+1] * refRatio(lev)[0];
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
