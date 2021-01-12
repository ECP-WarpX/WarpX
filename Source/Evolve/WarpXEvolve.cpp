/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Aurore Blelly
 *                     Axel Huebl, Burlen Loring, David Grote
 *                     Glenn Richardson, Jean-Luc Vay, Luca Fedeli
 *                     Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 *                     Weiqun Zhang, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "FieldSolver/WarpX_QED_K.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Python/WarpX_py.H"
#ifdef WARPX_USE_PSATD
#   include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif

#include <cmath>
#include <limits>


using namespace amrex;

void
WarpX::Evolve (int numsteps)
{
    WARPX_PROFILE("WarpX::Evolve()");

    Real cur_time = t_new[0];

    if (do_compute_max_step_from_zmax) {
        computeMaxStepBoostAccelerator(geom[0]);
    }

    int numsteps_max;
    if (numsteps < 0) {  // Note that the default argument is numsteps = -1
        numsteps_max = max_step;
    } else {
        numsteps_max = std::min(istep[0]+numsteps, max_step);
    }

    Real walltime, walltime_start = amrex::second();
    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {
        Real walltime_beg_step = amrex::second();

        multi_diags->NewIteration();

        // Start loop on time steps
        amrex::Print() << "\nSTEP " << step+1 << " starts ...\n";

        if (warpx_py_beforestep) warpx_py_beforestep();

        amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(0);
        if (cost) {
            if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD)
                amrex::Abort("LoadBalance for PSATD: TODO");
            if (step > 0 && load_balance_intervals.contains(step+1))
            {
                LoadBalance();

                // Reset the costs to 0
                ResetCosts();
            }
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                cost = WarpX::getCosts(lev);
                if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
                {
                    // Perform running average of the costs
                    // (Giving more importance to most recent costs; only needed
                    // for timers update, heuristic load balance considers the
                    // instantaneous costs)
                    for (int i : cost->IndexArray())
                    {
                        (*cost)[i] *= (1. - 2./load_balance_intervals.localPeriod(step+1));
                    }
                }
            }
        }

        // At the beginning, we have B^{n} and E^{n}.
        // Particles have p^{n} and x^{n}.
        // is_synchronized is true.
        if (is_synchronized) {
            // Not called at each iteration, so exchange all guard cells
            FillBoundaryE(guard_cells.ng_alloc_EB, guard_cells.ng_Extra);
            FillBoundaryB(guard_cells.ng_alloc_EB, guard_cells.ng_Extra);
            UpdateAuxilaryData();
            FillBoundaryAux(guard_cells.ng_UpdateAux);
            // on first step, push p by -0.5*dt
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                mypc->PushP(lev, -0.5_rt*dt[lev],
                            *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                            *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
            }
            is_synchronized = false;
        } else {
            // Beyond one step, we have E^{n} and B^{n}.
            // Particles have p^{n-1/2} and x^{n}.

            // E and B are up-to-date inside the domain only
            FillBoundaryE(guard_cells.ng_FieldGather, guard_cells.ng_Extra);
            FillBoundaryB(guard_cells.ng_FieldGather, guard_cells.ng_Extra);
            // E and B: enough guard cells to update Aux or call Field Gather in fp and cp
            // Need to update Aux on lower levels, to interpolate to higher levels.
            if (fft_do_time_averaging)
            {
                FillBoundaryE_avg(guard_cells.ng_FieldGather, guard_cells.ng_Extra);
                FillBoundaryB_avg(guard_cells.ng_FieldGather, guard_cells.ng_Extra);
            }
            // TODO Remove call to FillBoundaryAux before UpdateAuxilaryData?
            if (WarpX::maxwell_solver_id != MaxwellSolverAlgo::PSATD)
                FillBoundaryAux(guard_cells.ng_UpdateAux);
            UpdateAuxilaryData();
            FillBoundaryAux(guard_cells.ng_UpdateAux);
        }
        if (do_subcycling == 0 || finest_level == 0) {
            OneStep_nosub(cur_time);
            // E : guard cells are up-to-date
            // B : guard cells are NOT up-to-date
            // F : guard cells are NOT up-to-date
        } else if (do_subcycling == 1 && finest_level == 1) {
            OneStep_sub1(cur_time);
        } else {
            amrex::Print() << "Error: do_subcycling = " << do_subcycling << std::endl;
            amrex::Abort("Unsupported do_subcycling type");
        }

        if (num_mirrors>0){
            applyMirrors(cur_time);
            // E : guard cells are NOT up-to-date
            // B : guard cells are NOT up-to-date
        }

        if (warpx_py_beforeEsolve) warpx_py_beforeEsolve();

        if (cur_time + dt[0] >= stop_time - 1.e-3*dt[0] || step == numsteps_max-1) {
            // At the end of last step, push p by 0.5*dt to synchronize
            UpdateAuxilaryData();
            FillBoundaryAux(guard_cells.ng_UpdateAux);
            for (int lev = 0; lev <= finest_level; ++lev) {
                mypc->PushP(lev, 0.5_rt*dt[lev],
                            *Efield_aux[lev][0],*Efield_aux[lev][1],
                            *Efield_aux[lev][2],
                            *Bfield_aux[lev][0],*Bfield_aux[lev][1],
                            *Bfield_aux[lev][2]);
            }
            is_synchronized = true;
        }

        if (warpx_py_afterEsolve) warpx_py_afterEsolve();

        for (int lev = 0; lev <= max_level; ++lev) {
            ++istep[lev];
        }

        cur_time += dt[0];

        ShiftGalileanBoundary();

        if (do_back_transformed_diagnostics) {
            std::unique_ptr<MultiFab> cell_centered_data = nullptr;
            if (WarpX::do_back_transformed_fields) {
                cell_centered_data = GetCellCenteredData();
            }
            myBFD->writeLabFrameData(cell_centered_data.get(), *mypc, geom[0], cur_time, dt[0]);
        }

        bool move_j = is_synchronized;
        // If is_synchronized we need to shift j too so that next step we can evolve E by dt/2.
        // We might need to move j because we are going to make a plotfile.

        int num_moved = MoveWindow(move_j);

        mypc->ApplyBoundaryConditions();

        // Electrostatic solver: particles can move by an arbitrary number of cells
        if( do_electrostatic != ElectrostaticSolverAlgo::None )
        {
            mypc->Redistribute();
        } else
        {
            // Electromagnetic solver: due to CFL condition, particles can
            // only move by one or two cells per time step
            if (max_level == 0) {
                int num_redistribute_ghost = num_moved;
                if ((m_v_galilean[0]!=0) or (m_v_galilean[1]!=0) or (m_v_galilean[2]!=0)) {
                    // Galilean algorithm ; particles can move by up to 2 cells
                    num_redistribute_ghost += 2;
                } else {
                    // Standard algorithm ; particles can move by up to 1 cell
                    num_redistribute_ghost += 1;
                }
                mypc->RedistributeLocal(num_redistribute_ghost);
            }
            else {
                mypc->Redistribute();
            }
        }


        if (sort_intervals.contains(step+1)) {
            amrex::Print() << "re-sorting particles \n";
            mypc->SortParticlesByBin(sort_bin_size);
        }

        amrex::Print()<< "STEP " << step+1 << " ends." << " TIME = " << cur_time
                      << " DT = " << dt[0] << "\n";
        Real walltime_end_step = amrex::second();
        walltime = walltime_end_step - walltime_start;
        amrex::Print()<< "Walltime = " << walltime
                      << " s; This step = " << walltime_end_step-walltime_beg_step
                      << " s; Avg. per step = " << walltime/(step+1) << " s\n";

        // sync up time
        for (int i = 0; i <= max_level; ++i) {
            t_new[i] = cur_time;
        }

        /// reduced diags
        if (reduced_diags->m_plot_rd != 0)
        {
            reduced_diags->ComputeDiags(step);
            reduced_diags->WriteToFile(step);
        }
        multi_diags->FilterComputePackFlush( step );

        if (cur_time >= stop_time - 1.e-3*dt[0]) {
            break;
        }

        if (warpx_py_afterstep) warpx_py_afterstep();

        // End loop on time steps
    }

    multi_diags->FilterComputePackFlush( istep[0], true );

    if (do_back_transformed_diagnostics) {
        myBFD->Flush(geom[0]);
    }
}

/* /brief Perform one PIC iteration, without subcycling
*  i.e. all levels/patches use the same timestep (that of the finest level)
*  for the field advance and particle pusher.
*/
void
WarpX::OneStep_nosub (Real cur_time)
{

    if( do_electrostatic != ElectrostaticSolverAlgo::None ) {
        // Electrostatic solver:
        // For each species: deposit charge and add the associated space-charge
        // E and B field to the grid ; this is done at the beginning of the PIC
        // loop (i.e. immediately after a `Redistribute` and before particle
        // positions are pushed) so that the particles do not deposit out of bound
        bool const reset_fields = true;
        ComputeSpaceChargeField( reset_fields );
    }

    // Loop over species. For each ionizable species, create particles in
    // product species.
    doFieldIonization();

    mypc->doCollisions( cur_time );
#ifdef WARPX_QED
    mypc->doQEDSchwinger();
#endif
    // Push particle from x^{n} to x^{n+1}
    //               from p^{n-1/2} to p^{n+1/2}
    // Deposit current j^{n+1/2}
    // Deposit charge density rho^{n}
    if (warpx_py_particleinjection) warpx_py_particleinjection();
    if (warpx_py_particlescraper) warpx_py_particlescraper();
    if (warpx_py_beforedeposition) warpx_py_beforedeposition();
    PushParticlesandDepose(cur_time);

    if (warpx_py_afterdeposition) warpx_py_afterdeposition();

// TODO
// Apply current correction in Fourier space: for domain decomposition with local
// FFTs over guard cells, apply this before calling SyncCurrent
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        if (!fft_periodic_single_box && current_correction)
            amrex::Abort(
                    "\nCurrent correction does not guarantee charge conservation with local FFTs over guard cells:\n"
                    "set psatd.periodic_single_box_fft=1 too, in order to guarantee charge conservation");
        if (!fft_periodic_single_box && (WarpX::current_deposition_algo == CurrentDepositionAlgo::Vay))
            amrex::Abort(
                    "\nVay current deposition does not guarantee charge conservation with local FFTs over guard cells:\n"
                    "set psatd.periodic_single_box_fft=1 too, in order to guarantee charge conservation");
    }

#ifdef WARPX_QED
    doQEDEvents();
#endif

    // +1 is necessary here because value of step seen by user (first step is 1) is different than
    // value of step in code (first step is 0)
    mypc->doResampling(istep[0]+1);

    // Synchronize J and rho
    SyncCurrent();
    SyncRho();

    // Apply current correction in Fourier space: for periodic single-box global FFTs
    // without guard cells, apply this after calling SyncCurrent
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        if (fft_periodic_single_box && current_correction) CurrentCorrection();
        if (fft_periodic_single_box && (WarpX::current_deposition_algo == CurrentDepositionAlgo::Vay))
            VayDeposition();
    }


    // At this point, J is up-to-date inside the domain, and E and B are
    // up-to-date including enough guard cells for first step of the field
    // solve.

    // For extended PML: copy J from regular grid to PML, and damp J in PML
    if (do_pml && pml_has_particles) CopyJPML();
    if (do_pml && do_pml_j_damping) DampJPML();

    if( do_electrostatic == ElectrostaticSolverAlgo::None ) {
        // Electromagnetic solver:
        // Push E and B from {n} to {n+1}
        // (And update guard cells immediately afterwards)
        if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            if (use_hybrid_QED)
            {
                WarpX::Hybrid_QED_Push(dt);
                FillBoundaryE(guard_cells.ng_alloc_EB, guard_cells.ng_Extra);
            }
            PushPSATD(dt[0]);
            FillBoundaryE(guard_cells.ng_alloc_EB, guard_cells.ng_Extra);
            FillBoundaryB(guard_cells.ng_alloc_EB, guard_cells.ng_Extra);

            if (use_hybrid_QED)
            {
                WarpX::Hybrid_QED_Push(dt);
                FillBoundaryE(guard_cells.ng_alloc_EB, guard_cells.ng_Extra);
            }
            if (do_pml) DampPML();
        } else {
            EvolveF(0.5_rt * dt[0], DtType::FirstHalf);
            FillBoundaryF(guard_cells.ng_FieldSolverF);
            EvolveB(0.5_rt * dt[0]); // We now have B^{n+1/2}

            FillBoundaryB(guard_cells.ng_FieldSolver, IntVect::TheZeroVector());
            if (WarpX::em_solver_medium == MediumForEM::Vacuum) {
                // vacuum medium
                EvolveE(dt[0]); // We now have E^{n+1}
            } else if (WarpX::em_solver_medium == MediumForEM::Macroscopic) {
                // macroscopic medium
                MacroscopicEvolveE(dt[0]); // We now have E^{n+1}
            } else {
                amrex::Abort(" Medium for EM is unknown \n");
            }

            FillBoundaryE(guard_cells.ng_FieldSolver, IntVect::TheZeroVector());
            EvolveF(0.5_rt * dt[0], DtType::SecondHalf);
            EvolveB(0.5_rt * dt[0]); // We now have B^{n+1}
            if (do_pml) {
                FillBoundaryF(guard_cells.ng_alloc_F);
                DampPML();
                FillBoundaryE(guard_cells.ng_MovingWindow, IntVect::TheZeroVector());
                FillBoundaryF(guard_cells.ng_MovingWindow);
                FillBoundaryB(guard_cells.ng_MovingWindow, IntVect::TheZeroVector());
            }
            // E and B are up-to-date in the domain, but all guard cells are
            // outdated.
            if (safe_guard_cells)
                FillBoundaryB(guard_cells.ng_alloc_EB, guard_cells.ng_Extra);
        } // !PSATD
    } // !do_electrostatic
}

/* /brief Perform one PIC iteration, with subcycling
*  i.e. The fine patch uses a smaller timestep (and steps more often)
*  than the coarse patch, for the field advance and particle pusher.
*
* This version of subcycling only works for 2 levels and with a refinement
* ratio of 2.
* The particles and fields of the fine patch are pushed twice
* (with dt[coarse]/2) in this routine.
* The particles of the coarse patch and mother grid are pushed only once
* (with dt[coarse]). The fields on the coarse patch and mother grid
* are pushed in a way which is equivalent to pushing once only, with
* a current which is the average of the coarse + fine current at the 2
* steps of the fine grid.
*
*/


void
WarpX::OneStep_sub1 (Real curtime)
{
    if( do_electrostatic != ElectrostaticSolverAlgo::None )
    {
        amrex::Abort("Electrostatic solver cannot be used with sub-cycling.");
    }

    // TODO: we could save some charge depositions
    // Loop over species. For each ionizable species, create particles in
    // product species.
    doFieldIonization();

#ifdef WARPX_QED
    doQEDEvents();
#endif

    // +1 is necessary here because value of step seen by user (first step is 1) is different than
    // value of step in code (first step is 0)
    mypc->doResampling(istep[0]+1);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(finest_level == 1, "Must have exactly two levels");
    const int fine_lev = 1;
    const int coarse_lev = 0;

    // i) Push particles and fields on the fine patch (first fine step)
    PushParticlesandDepose(fine_lev, curtime, DtType::FirstHalf);
    RestrictCurrentFromFineToCoarsePatch(fine_lev);
    RestrictRhoFromFineToCoarsePatch(fine_lev);
    ApplyFilterandSumBoundaryJ(fine_lev, PatchType::fine);
    NodalSyncJ(fine_lev, PatchType::fine);
    ApplyFilterandSumBoundaryRho(fine_lev, PatchType::fine, 0, 2*ncomps);
    NodalSyncRho(fine_lev, PatchType::fine, 0, 2);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);
    FillBoundaryF(fine_lev, PatchType::fine, guard_cells.ng_alloc_F);

    EvolveE(fine_lev, PatchType::fine, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::fine, guard_cells.ng_FieldGather);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::SecondHalf);

    if (do_pml) {
        FillBoundaryF(fine_lev, PatchType::fine, guard_cells.ng_alloc_F);
        DampPML(fine_lev, PatchType::fine);
        FillBoundaryE(fine_lev, PatchType::fine, guard_cells.ng_FieldGather);
    }

    FillBoundaryB(fine_lev, PatchType::fine, guard_cells.ng_FieldGather);

    // ii) Push particles on the coarse patch and mother grid.
    // Push the fields on the coarse patch and mother grid
    // by only half a coarse step (first half)
    PushParticlesandDepose(coarse_lev, curtime, DtType::Full);
    StoreCurrent(coarse_lev);
    AddCurrentFromFineLevelandSumBoundary(coarse_lev);
    AddRhoFromFineLevelandSumBoundary(coarse_lev, 0, ncomps);

    EvolveB(fine_lev, PatchType::coarse, dt[fine_lev]);
    EvolveF(fine_lev, PatchType::coarse, dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::coarse, guard_cells.ng_FieldGather);
    FillBoundaryF(fine_lev, PatchType::coarse, guard_cells.ng_FieldSolverF);

    EvolveE(fine_lev, PatchType::coarse, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::coarse, guard_cells.ng_FieldGather);

    EvolveB(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev]);
    EvolveF(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev], DtType::FirstHalf);
    FillBoundaryB(coarse_lev, PatchType::fine, guard_cells.ng_FieldGather + guard_cells.ng_Extra);
    FillBoundaryF(coarse_lev, PatchType::fine, guard_cells.ng_FieldSolverF);

    EvolveE(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev]);
    FillBoundaryE(coarse_lev, PatchType::fine, guard_cells.ng_FieldGather + guard_cells.ng_Extra);

    // TODO Remove call to FillBoundaryAux before UpdateAuxilaryData?
    FillBoundaryAux(guard_cells.ng_UpdateAux);
    // iii) Get auxiliary fields on the fine grid, at dt[fine_lev]
    UpdateAuxilaryData();
    FillBoundaryAux(guard_cells.ng_UpdateAux);

    // iv) Push particles and fields on the fine patch (second fine step)
    PushParticlesandDepose(fine_lev, curtime+dt[fine_lev], DtType::SecondHalf);
    RestrictCurrentFromFineToCoarsePatch(fine_lev);
    RestrictRhoFromFineToCoarsePatch(fine_lev);
    ApplyFilterandSumBoundaryJ(fine_lev, PatchType::fine);
    NodalSyncJ(fine_lev, PatchType::fine);
    ApplyFilterandSumBoundaryRho(fine_lev, PatchType::fine, 0, ncomps);
    NodalSyncRho(fine_lev, PatchType::fine, 0, 2);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);
    FillBoundaryF(fine_lev, PatchType::fine, guard_cells.ng_FieldSolverF);

    EvolveE(fine_lev, PatchType::fine, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::SecondHalf);

    if (do_pml) {
        DampPML(fine_lev, PatchType::fine);
        FillBoundaryE(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);
    }

    if ( safe_guard_cells )
        FillBoundaryF(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);
    FillBoundaryB(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);

    // v) Push the fields on the coarse patch and mother grid
    // by only half a coarse step (second half)
    RestoreCurrent(coarse_lev);
    AddCurrentFromFineLevelandSumBoundary(coarse_lev);
    AddRhoFromFineLevelandSumBoundary(coarse_lev, ncomps, ncomps);

    EvolveE(fine_lev, PatchType::coarse, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::coarse, guard_cells.ng_FieldSolver);

    EvolveB(fine_lev, PatchType::coarse, dt[fine_lev]);
    EvolveF(fine_lev, PatchType::coarse, dt[fine_lev], DtType::SecondHalf);

    if (do_pml) {
        FillBoundaryF(fine_lev, PatchType::fine, guard_cells.ng_FieldSolverF);
        DampPML(fine_lev, PatchType::coarse); // do it twice
        DampPML(fine_lev, PatchType::coarse);
        FillBoundaryE(fine_lev, PatchType::coarse, guard_cells.ng_alloc_EB);
    }

    FillBoundaryB(fine_lev, PatchType::coarse, guard_cells.ng_FieldSolver);

    FillBoundaryF(fine_lev, PatchType::coarse, guard_cells.ng_FieldSolverF);

    EvolveE(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev]);
    FillBoundaryE(coarse_lev, PatchType::fine, guard_cells.ng_FieldSolver);

    EvolveB(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev]);
    EvolveF(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev], DtType::SecondHalf);

    if (do_pml) {
        if (do_moving_window){
            // Exchance guard cells of PMLs only (0 cells are exchanged for the
            // regular B field MultiFab). This is required as B and F have just been
            // evolved.
            FillBoundaryB(coarse_lev, PatchType::fine, IntVect::TheZeroVector());
            FillBoundaryF(coarse_lev, PatchType::fine, IntVect::TheZeroVector());
        }
        DampPML(coarse_lev, PatchType::fine);
        if ( safe_guard_cells )
            FillBoundaryE(coarse_lev, PatchType::fine, guard_cells.ng_FieldSolver);
    }
    if ( safe_guard_cells )
        FillBoundaryB(coarse_lev, PatchType::fine, guard_cells.ng_FieldSolver);
}

void
WarpX::doFieldIonization ()
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        doFieldIonization(lev);
    }
}

void
WarpX::doFieldIonization (int lev)
{
    mypc->doFieldIonization(lev,
                            *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                            *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
}

#ifdef WARPX_QED
void
WarpX::doQEDEvents ()
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        doQEDEvents(lev);
    }
}

void
WarpX::doQEDEvents (int lev)
{
    mypc->doQedEvents(lev,
                      *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                      *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
}
#endif

void
WarpX::PushParticlesandDepose (amrex::Real cur_time)
{
    // Evolve particles to p^{n+1/2} and x^{n+1}
    // Depose current, j^{n+1/2}
    for (int lev = 0; lev <= finest_level; ++lev) {
        PushParticlesandDepose(lev, cur_time);
    }
}

void
WarpX::PushParticlesandDepose (int lev, amrex::Real cur_time, DtType a_dt_type)
{
    mypc->Evolve(lev,
                 *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                 *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2],
                 *Efield_avg_aux[lev][0],*Efield_avg_aux[lev][1],*Efield_avg_aux[lev][2],
                 *Bfield_avg_aux[lev][0],*Bfield_avg_aux[lev][1],*Bfield_avg_aux[lev][2],
                 *current_fp[lev][0],*current_fp[lev][1],*current_fp[lev][2],
                 current_buf[lev][0].get(), current_buf[lev][1].get(), current_buf[lev][2].get(),
                 rho_fp[lev].get(), charge_buf[lev].get(),
                 Efield_cax[lev][0].get(), Efield_cax[lev][1].get(), Efield_cax[lev][2].get(),
                 Bfield_cax[lev][0].get(), Bfield_cax[lev][1].get(), Bfield_cax[lev][2].get(),
                 cur_time, dt[lev], a_dt_type);
#ifdef WARPX_DIM_RZ
    // This is called after all particles have deposited their current and charge.
    ApplyInverseVolumeScalingToCurrentDensity(current_fp[lev][0].get(), current_fp[lev][1].get(), current_fp[lev][2].get(), lev);
    if (current_buf[lev][0].get()) {
        ApplyInverseVolumeScalingToCurrentDensity(current_buf[lev][0].get(), current_buf[lev][1].get(), current_buf[lev][2].get(), lev-1);
    }
    if (rho_fp[lev].get()) {
        ApplyInverseVolumeScalingToChargeDensity(rho_fp[lev].get(), lev);
        if (charge_buf[lev].get()) {
            ApplyInverseVolumeScalingToChargeDensity(charge_buf[lev].get(), lev-1);
        }
    }
#endif
}

/* \brief computes max_step for wakefield simulation in boosted frame.
 * \param geom: Geometry object that contains simulation domain.
 *
 * max_step is set so that the simulation stop when the lower corner of the
 * simulation box passes input parameter zmax_plasma_to_compute_max_step.
 */
void
WarpX::computeMaxStepBoostAccelerator(const amrex::Geometry& a_geom){
    // Sanity checks: can use zmax_plasma_to_compute_max_step only if
    // the moving window and the boost are all in z direction.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::moving_window_dir == AMREX_SPACEDIM-1,
        "Can use zmax_plasma_to_compute_max_step only if " +
        "moving window along z. TODO: all directions.");
    if (gamma_boost > 1){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            (WarpX::boost_direction[0]-0)*(WarpX::boost_direction[0]-0) +
            (WarpX::boost_direction[1]-0)*(WarpX::boost_direction[1]-0) +
            (WarpX::boost_direction[2]-1)*(WarpX::boost_direction[2]-1) < 1.e-12,
            "Can use zmax_plasma_to_compute_max_step in boosted frame only if " +
            "warpx.boost_direction = z. TODO: all directions.");
    }

    // Lower end of the simulation domain. All quantities are given in boosted
    // frame except zmax_plasma_to_compute_max_step.
    const Real zmin_domain_boost = a_geom.ProbLo(AMREX_SPACEDIM-1);
    // End of the plasma: Transform input argument
    // zmax_plasma_to_compute_max_step to boosted frame.
    const Real len_plasma_boost = zmax_plasma_to_compute_max_step/gamma_boost;
    // Plasma velocity
    const Real v_plasma_boost = -beta_boost * PhysConst::c;
    // Get time at which the lower end of the simulation domain passes the
    // upper end of the plasma (in the z direction).
    const Real interaction_time_boost = (len_plasma_boost-zmin_domain_boost)/
        (moving_window_v-v_plasma_boost);
    // Divide by dt, and update value of max_step.
    int computed_max_step;
    if (do_subcycling){
        computed_max_step = static_cast<int>(interaction_time_boost/dt[0]);
    } else {
        computed_max_step =
            static_cast<int>(interaction_time_boost/dt[maxLevel()]);
    }
    max_step = computed_max_step;
    Print()<<"max_step computed in computeMaxStepBoostAccelerator: "
           <<computed_max_step<<std::endl;
}

/* \brief Apply perfect mirror condition inside the box (not at a boundary).
 * In practice, set all fields to 0 on a section of the simulation domain
 * (as for a perfect conductor with a given thickness).
 * The mirror normal direction has to be parallel to the z axis.
 */
void
WarpX::applyMirrors(Real time){
    // Loop over the mirrors
    for(int i_mirror=0; i_mirror<num_mirrors; ++i_mirror){
        // Get mirror properties (lower and upper z bounds)
        Real z_min = mirror_z[i_mirror];
        Real z_max_tmp = z_min + mirror_z_width[i_mirror];
        // Boost quantities for boosted frame simulations
        if (gamma_boost>1){
            z_min = z_min/gamma_boost - PhysConst::c*beta_boost*time;
            z_max_tmp = z_max_tmp/gamma_boost - PhysConst::c*beta_boost*time;
        }
        // Loop over levels
        for(int lev=0; lev<=finest_level; lev++){
            // Make sure that the mirror contains at least
            // mirror_z_npoints[i_mirror] cells
            Real dz = WarpX::CellSize(lev)[2];
            Real z_max = std::max(z_max_tmp,
                                  z_min+mirror_z_npoints[i_mirror]*dz);
            // Get fine patch field MultiFabs
            MultiFab& Ex = *Efield_fp[lev][0].get();
            MultiFab& Ey = *Efield_fp[lev][1].get();
            MultiFab& Ez = *Efield_fp[lev][2].get();
            MultiFab& Bx = *Bfield_fp[lev][0].get();
            MultiFab& By = *Bfield_fp[lev][1].get();
            MultiFab& Bz = *Bfield_fp[lev][2].get();
            // Set each field to zero between z_min and z_max
            NullifyMF(Ex, lev, z_min, z_max);
            NullifyMF(Ey, lev, z_min, z_max);
            NullifyMF(Ez, lev, z_min, z_max);
            NullifyMF(Bx, lev, z_min, z_max);
            NullifyMF(By, lev, z_min, z_max);
            NullifyMF(Bz, lev, z_min, z_max);
            if (lev>0){
                // Get coarse patch field MultiFabs
                MultiFab& cEx = *Efield_cp[lev][0].get();
                MultiFab& cEy = *Efield_cp[lev][1].get();
                MultiFab& cEz = *Efield_cp[lev][2].get();
                MultiFab& cBx = *Bfield_cp[lev][0].get();
                MultiFab& cBy = *Bfield_cp[lev][1].get();
                MultiFab& cBz = *Bfield_cp[lev][2].get();
                // Set each field to zero between z_min and z_max
                NullifyMF(cEx, lev, z_min, z_max);
                NullifyMF(cEy, lev, z_min, z_max);
                NullifyMF(cEz, lev, z_min, z_max);
                NullifyMF(cBx, lev, z_min, z_max);
                NullifyMF(cBy, lev, z_min, z_max);
                NullifyMF(cBz, lev, z_min, z_max);
            }
        }
    }
}

// Apply current correction in Fourier space
void
WarpX::CurrentCorrection ()
{
#ifdef WARPX_USE_PSATD
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD)
    {
        for ( int lev = 0; lev <= finest_level; ++lev )
        {
            spectral_solver_fp[lev]->CurrentCorrection( current_fp[lev], rho_fp[lev] );
            if ( spectral_solver_cp[lev] ) spectral_solver_cp[lev]->CurrentCorrection( current_cp[lev], rho_cp[lev] );
        }
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
            "WarpX::CurrentCorrection: only implemented for spectral solver.");
    }
#else
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
    "WarpX::CurrentCorrection: requires WarpX build with spectral solver support.");
#endif
}

// Compute current from Vay deposition in Fourier space
void
WarpX::VayDeposition ()
{
#ifdef WARPX_USE_PSATD
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            spectral_solver_fp[lev]->VayDeposition(current_fp[lev]);
            if (spectral_solver_cp[lev]) spectral_solver_cp[lev]->VayDeposition(current_cp[lev]);
        }
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
            "WarpX::VayDeposition: only implemented for spectral solver.");
    }
#else
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
    "WarpX::CurrentCorrection: requires WarpX build with spectral solver support.");
#endif
}
