/* Copyright 2020
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Utils/WarpXAlgorithmSelection.H"
#ifndef WARPX_DIM_RZ
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#endif

/**
 * Determine the timestep of the simulation. */
void
WarpX::ComputeDt ()
{
    // Determine
    const amrex::Real* dx = geom[max_level].CellSize();
    amrex::Real deltat = 0.;

#if WARPX_USE_PSATD
    // Computation of dt for spectral algorithm

#   if (defined WARPX_DIM_RZ)
    // - In RZ geometry: dz/c
    deltat = cfl * dx[1]/PhysConst::c;
#   elif (defined WARPX_DIM_XZ)
    // - In Cartesian 2D geometry: determined by the minimum cell size in all direction
    deltat = cfl * std::min( dx[0], dx[1] )/PhysConst::c;
#   else
    // - In Cartesian 3D geometry: determined by the minimum cell size in all direction
    deltat = cfl * std::min( dx[0], std::min( dx[1], dx[2] ) )/PhysConst::c;
#   endif


#else
    // Computation of dt for FDTD algorithm

#   ifdef WARPX_DIM_RZ
    // In the rz case, the Courant limit has been evaluated
    // semi-analytically by R. Lehe, and resulted in the following
    // coefficients.
    // NB : Here the coefficient for m=1 as compared to this document,
    // as it was observed in practice that this coefficient was not
    // high enough (The simulation became unstable).
    amrex::Real multimode_coeffs[6] = { 0.2105, 1.0, 3.5234, 8.5104, 15.5059, 24.5037 };
    amrex::Real multimode_alpha;
    if (n_rz_azimuthal_modes < 7) {
        // Use the table of the coefficients
        multimode_alpha = multimode_coeffs[n_rz_azimuthal_modes-1];
    } else {
        // Use a realistic extrapolation
        multimode_alpha = (n_rz_azimuthal_modes - 1)*(n_rz_azimuthal_modes - 1) - 0.4;
    }
    deltat  = cfl * 1./( std::sqrt((1+multimode_alpha)/(dx[0]*dx[0]) + 1./(dx[1]*dx[1])) * PhysConst::c );
#   else
    // - In Cartesian geometry
    if (do_nodal) {
        deltat = cfl * CartesianNodalAlgorithm::ComputeMaxDt( dx );
    } else if (maxwell_fdtd_solver_id == MaxwellSolverAlgo::Yee) {
        deltat = cfl * CartesianYeeAlgorithm::ComputeMaxDt( dx );
    } else if (maxwell_fdtd_solver_id == MaxwellSolverAlgo::CKC) {
        deltat = cfl * CartesianCKCAlgorithm::ComputeMaxDt( dx );
    } else {
        amrex::Abort("Unknown algorithm");
    }
#   endif

#endif

    dt.resize(0);
    dt.resize(max_level+1,deltat);

    if (do_subcycling) {
        for (int lev = max_level-1; lev >= 0; --lev) {
            dt[lev] = dt[lev+1] * refRatio(lev)[0];
        }
    }

    if (do_electrostatic) {
        dt[0] = const_dt;
    }

    for (int lev=0; lev <= max_level; lev++) {
        const amrex::Real* dx_lev = geom[lev].CellSize();
        amrex::Print()<<"Level "<<lev<<": dt = "<<dt[lev]
               <<" ; dx = "<<dx_lev[0]
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
               <<" ; dz = "<<dx_lev[1]<<'\n';
#elif (defined WARPX_DIM_3D)
               <<" ; dy = "<<dx_lev[1]
               <<" ; dz = "<<dx_lev[2]<<'\n';
#endif
    }
}
