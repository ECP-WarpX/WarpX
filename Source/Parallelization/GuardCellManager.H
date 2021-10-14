/* Copyright 2019-2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef GUARDCELLMANAGER_H_
#define GUARDCELLMANAGER_H_

#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

/**
 * \brief This class computes and stores the number of guard cells needed for
 * the allocation of the MultiFabs and required for each part of the PIC loop.
 */
class guardCellManager{

public:

    /**
     * \brief Initialize number of guard cells depending on the options used.
     *
     * \param do_subcycling bool, whether to use subcycling
     * \param do_fdtd_nci_corr bool, whether to use Godfrey NCI corrector
     * \param do_nodal bool, whether the field solver is nodal
     * \param do_moving_window bool, whether to use moving window
     * \param aux_is_nodal bool, true if the aux grid is nodal
     * \param moving_window_dir direction of moving window
     * \param nox order of current deposition
     * \param nox_fft order of PSATD in x direction
     * \param noy_fft order of PSATD in y direction
     * \param noz_fft order of PSATD in z direction
     * \param nci_corr_stencil stencil of NCI corrector
     * \param maxwell_solver_id if of Maxwell solver
     * \param max_level max level of the simulation
     */
    void Init(
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
        const amrex::Vector<amrex::Real> v_galilean,
        const amrex::Vector<amrex::Real> v_comoving,
        const bool safe_guard_cells,
        const int do_electrostatic,
        const int do_multi_J,
        const bool fft_do_time_averaging);

    // Guard cells allocated for MultiFabs E and B
    amrex::IntVect ng_alloc_EB = amrex::IntVect::TheZeroVector();
    // Guard cells allocated for MultiFab J
    amrex::IntVect ng_alloc_J = amrex::IntVect::TheZeroVector();
    // Guard cells allocated for MultiFab Rho
    amrex::IntVect ng_alloc_Rho = amrex::IntVect::TheZeroVector();
    // Guard cells allocated for MultiFab F
    amrex::IntVect ng_alloc_F = amrex::IntVect::TheZeroVector();
    // Guard cells allocated for MultiFab G
    amrex::IntVect ng_alloc_G = amrex::IntVect::TheZeroVector();

    // Guard cells exchanged for specific parts of the PIC loop

    // Number of guard cells of E and B that must exchanged before Field Solver
    amrex::IntVect ng_FieldSolver = amrex::IntVect::TheZeroVector();
    // Number of guard cells of F that must exchanged before Field Solver
    amrex::IntVect ng_FieldSolverF = amrex::IntVect::TheZeroVector();
    // Number of guard cells of G that must be exchanged before Field Solver
    amrex::IntVect ng_FieldSolverG = amrex::IntVect::TheZeroVector();
    // Number of guard cells of E and B that must exchanged before Field Gather
    amrex::IntVect ng_FieldGather = amrex::IntVect::TheZeroVector();
    // Number of guard cells of E and B that must exchanged before updating the Aux grid
    amrex::IntVect ng_UpdateAux = amrex::IntVect::TheZeroVector();
    // Number of guard cells of all MultiFabs that must exchanged before moving window
    amrex::IntVect ng_MovingWindow = amrex::IntVect::TheZeroVector();
    // Number of guard cells of E and B that are exchanged immediatly after the main PSATD push
    amrex::IntVect ng_afterPushPSATD = amrex::IntVect::TheZeroVector();

    // Number of guard cells for local deposition of J and rho
    amrex::IntVect ng_depos_J   = amrex::IntVect::TheZeroVector();
    amrex::IntVect ng_depos_rho = amrex::IntVect::TheZeroVector();
};

#endif // GUARDCELLMANAGER_H_
