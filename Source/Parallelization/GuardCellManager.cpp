/* Copyright 2019-2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "GuardCellManager.H"
#include "Filter/NCIGodfreyFilter.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"

#include <AMReX_ParmParse.H>
#include <AMReX.H>

using namespace amrex;

void
guardCellManager::Init (
    const amrex::Real dt,
    const amrex::RealVect dx,
    const bool do_subcycling,
    const bool do_fdtd_nci_corr,
    const bool do_nodal,
    const bool do_moving_window,
    const int moving_window_dir,
    const int nox,
    const int nox_fft, const int noy_fft, const int noz_fft,
    const int nci_corr_stencil,
    const int maxwell_solver_id,
    const int max_level,
    const amrex::Array<amrex::Real,3> v_galilean,
    const amrex::Array<amrex::Real,3> v_comoving,
    const bool safe_guard_cells,
    const int do_electrostatic)
{
    // When using subcycling, the particles on the fine level perform two pushes
    // before being redistributed ; therefore, we need one extra guard cell
    // (the particles may move by 2*c*dt)
    int ngx_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;
    int ngy_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;
    int ngz_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;

    const bool galilean = (v_galilean[0] != 0. || v_galilean[1] != 0. || v_galilean[2] != 0.);
    const bool comoving = (v_comoving[0] != 0. || v_comoving[1] != 0. || v_comoving[2] != 0.);

    // Add one guard cell in the case of the Galilean or comoving algorithms
    if (galilean || comoving)
    {
      ngx_tmp += 1;
      ngy_tmp += 1;
      ngz_tmp += 1;
    }

    // Ex, Ey, Ez, Bx, By, and Bz have the same number of ghost cells.
    // jx, jy, jz and rho have the same number of ghost cells.
    // E and B have the same number of ghost cells as j and rho if NCI filter is not used,
    // but different number of ghost cells in z-direction if NCI filter is used.
    // The number of cells should be even, in order to easily perform the
    // interpolation from coarse grid to fine grid.
    int ngx = (ngx_tmp % 2) ? ngx_tmp+1 : ngx_tmp;  // Always even number
    int ngy = (ngy_tmp % 2) ? ngy_tmp+1 : ngy_tmp;  // Always even number
    int ngz_nonci = (ngz_tmp % 2) ? ngz_tmp+1 : ngz_tmp;  // Always even number
    int ngz;
    if (do_fdtd_nci_corr) {
        int ng = ngz_tmp + nci_corr_stencil;
        ngz = (ng % 2) ? ng+1 : ng;
    } else {
        ngz = ngz_nonci;
    }

    // J is only interpolated from fine to coarse (not coarse to fine)
    // and therefore does not need to be even.
    int ngJx = ngx_tmp;
    int ngJy = ngy_tmp;
    int ngJz = ngz_tmp;

    // When calling the moving window (with one level of refinement),  we shift
    // the fine grid by 2 cells ; therefore, we need at least 2 guard cells
    // on level 1. This may not be necessary for level 0.
    if (do_moving_window) {
        ngx = std::max(ngx,2);
        ngy = std::max(ngy,2);
        ngz = std::max(ngz,2);
        ngJx = std::max(ngJx,2);
        ngJy = std::max(ngJy,2);
        ngJz = std::max(ngJz,2);
    }

#if (AMREX_SPACEDIM == 3)
    ng_alloc_EB = IntVect(ngx,ngy,ngz);
    ng_alloc_J = IntVect(ngJx,ngJy,ngJz);
#elif (AMREX_SPACEDIM == 2)
    ng_alloc_EB = IntVect(ngx,ngz);
    ng_alloc_J = IntVect(ngJx,ngJz);
#endif

    // TODO Adding one cell for rho should not be necessary, given that the number of guard cells
    // now takes into account the time step (see code block below). However, this does seem to be
    // necessary in order to avoid some remaining instances of out-of-bound array access in
    // simulations with large time steps (revealed by building WarpX with BOUND_CHECK = TRUE).
    ng_alloc_Rho = ng_alloc_J+1;

    // Electromagnetic simulations: account for change in particle positions within half a time step
    // for current deposition and within one time step for charge deposition (since rho is needed
    // both at the beginning and at the end of the PIC iteration)
    if (do_electrostatic == ElectrostaticSolverAlgo::None)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
            ng_alloc_Rho[i] += static_cast<int>(std::ceil(PhysConst::c * dt / dx[i]));
            ng_alloc_J[i]   += static_cast<int>(std::ceil(PhysConst::c * dt / dx[i] * 0.5_rt));
        }
    }

    // Number of guard cells for local deposition of J and rho
    ng_depos_J   = ng_alloc_J;
    ng_depos_rho = ng_alloc_Rho;

    // After pushing particle
    int ng_alloc_F_int = (do_moving_window) ? 2 : 0;
    // CKC solver requires one additional guard cell
    if (maxwell_solver_id == MaxwellSolverAlgo::CKC) ng_alloc_F_int = std::max( ng_alloc_F_int, 1 );
    ng_alloc_F = IntVect(AMREX_D_DECL(ng_alloc_F_int, ng_alloc_F_int, ng_alloc_F_int));

    if (maxwell_solver_id == MaxwellSolverAlgo::PSATD)
    {
        // The number of guard cells should be enough to contain the stencil of the FFT solver.
        //
        // Here, this number (ngFFT) is determined empirically to be the order of the solver or
        // half the order of the solver, depending on other various numerical parameters.
        //
        // With the standard PSATD algorithm, simulations on staggered grids usually work fine
        // with a number of guard cells equal to half the number of guard cells that would be
        // used on nodal grids, in all directions x, y and z.
        //
        // On the other hand, with the Galilean PSATD or averaged Galilean PSATD algorithms,
        // with a Galilean coordinate transformation directed only in z, it seems more robust
        // to set the same number of guard cells in z, irrespective of whether the simulation
        // runs on nodal grids or staggered grids (typically with centering of fields and/or
        // currents in the latter case). This does not seem to be necessary in x and y,
        // where it still seems fine to set half the number of guard cells of the nodal case.

        int ngFFt_x = do_nodal ? nox_fft : nox_fft / 2;
        int ngFFt_y = do_nodal ? noy_fft : noy_fft / 2;
        int ngFFt_z = (do_nodal || galilean) ? noz_fft : noz_fft / 2;

        ParmParse pp_psatd("psatd");
        pp_psatd.query("nx_guard", ngFFt_x);
        pp_psatd.query("ny_guard", ngFFt_y);
        pp_psatd.query("nz_guard", ngFFt_z);

#if (AMREX_SPACEDIM == 3)
        IntVect ngFFT = IntVect(ngFFt_x, ngFFt_y, ngFFt_z);
#elif (AMREX_SPACEDIM == 2)
        IntVect ngFFT = IntVect(ngFFt_x, ngFFt_z);
#endif

        // All boxes should have the same number of guard cells, to avoid temporary parallel copies:
        // thus we take the maximum of the required number of guard cells over all available fields.
        for (int i_dim = 0; i_dim < AMREX_SPACEDIM; i_dim++) {
            int ng_required = ngFFT[i_dim];
            // Get the max
            ng_required = std::max(ng_required, ng_alloc_EB[i_dim]);
            ng_required = std::max(ng_required, ng_alloc_J[i_dim]);
            ng_required = std::max(ng_required, ng_alloc_Rho[i_dim]);
            ng_required = std::max(ng_required, ng_alloc_F[i_dim]);
            // Set the guard cells to this max
            ng_alloc_EB[i_dim] = ng_required;
            ng_alloc_J[i_dim] = ng_required;
            ng_alloc_F[i_dim] = ng_required;
            ng_alloc_Rho[i_dim] = ng_required;
            ng_alloc_F_int = ng_required;
        }
        ng_alloc_F = IntVect(AMREX_D_DECL(ng_alloc_F_int, ng_alloc_F_int, ng_alloc_F_int));
    }

    // Compute number of cells required for Field Solver
    if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        ng_FieldSolver = ng_alloc_EB;
        ng_FieldSolverF = ng_alloc_EB;
    } else {
        ng_FieldSolver = IntVect(AMREX_D_DECL(1, 1, 1));
        ng_FieldSolverF = IntVect(AMREX_D_DECL(1, 1, 1));
    }

    if (safe_guard_cells){
        // Run in safe mode: exchange all allocated guard cells at each
        // call of FillBoundary
        ng_FieldSolver = ng_alloc_EB;
        ng_FieldSolverF = ng_alloc_F;
        ng_FieldGather = ng_alloc_EB;
        ng_UpdateAux = ng_alloc_EB;
        if (do_moving_window){
            ng_MovingWindow = ng_alloc_EB;
        }
    } else {

        ng_FieldSolver = ng_FieldSolver.min(ng_alloc_EB);

        // Compute number of cells required for Field Gather
        int FGcell[4] = {0,1,1,2}; // Index is nox
        IntVect ng_FieldGather_noNCI = IntVect(AMREX_D_DECL(FGcell[nox],FGcell[nox],FGcell[nox]));
        ng_FieldGather_noNCI = ng_FieldGather_noNCI.min(ng_alloc_EB);
        // If NCI filter, add guard cells in the z direction
        IntVect ng_NCIFilter = IntVect::TheZeroVector();
        if (do_fdtd_nci_corr)
            ng_NCIFilter[AMREX_SPACEDIM-1] = NCIGodfreyFilter::m_stencil_width;
        // Note: communications of guard cells for bilinear filter are handled
        // separately.
        ng_FieldGather = ng_FieldGather_noNCI + ng_NCIFilter;

        // Guard cells for auxiliary grid.
        // Not sure why there is a 2* here...
        ng_UpdateAux = 2*ng_FieldGather_noNCI + ng_NCIFilter;

        // Make sure we do not exchange more guard cells than allocated.
        ng_FieldGather = ng_FieldGather.min(ng_alloc_EB);
        ng_UpdateAux = ng_UpdateAux.min(ng_alloc_EB);
        ng_FieldSolverF = ng_FieldSolverF.min(ng_alloc_F);
        // Only FillBoundary(ng_FieldGather) is called between consecutive
        // field solves. So ng_FieldGather must have enough cells
        // for the field solve too.
        ng_FieldGather = ng_FieldGather.max(ng_FieldSolver);

        if (do_moving_window){
            ng_MovingWindow[moving_window_dir] = 1;
        }
    }
}
