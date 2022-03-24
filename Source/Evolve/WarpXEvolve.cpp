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

#include "BoundaryConditions/PML.H"
#include "Diagnostics/BackTransformedDiagnostic.H"
#include "Diagnostics/MultiDiagnostics.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "Evolve/WarpXDtType.H"
#ifdef WARPX_USE_PSATD
#   ifdef WARPX_DIM_RZ
#       include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#   else
#       include "FieldSolver/SpectralSolver/SpectralSolver.H"
#   endif
#endif
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Python/WarpX_py.H"
#include "Utils/IntervalsParser.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <vector>

using namespace amrex;

void
WarpX::Evolve (int numsteps)
{
    WARPX_PROFILE_REGION("WarpX::Evolve()");
    WARPX_PROFILE("WarpX::Evolve()");

    Real cur_time = t_new[0];

    int numsteps_max;
    if (numsteps < 0) {  // Note that the default argument is numsteps = -1
        numsteps_max = max_step;
    } else {
        numsteps_max = std::min(istep[0]+numsteps, max_step);
    }

    bool early_params_checked = false; // check typos in inputs after step 1

    static Real evolve_time = 0;

    const int step_begin = istep[0];
    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {
        WARPX_PROFILE("WarpX::Evolve::step");
        Real evolve_time_beg_step = amrex::second();

        multi_diags->NewIteration();

        // Start loop on time steps
        if (verbose) {
            amrex::Print() << "\nSTEP " << step+1 << " starts ...\n";
        }
        ExecutePythonCallback("beforestep");

        amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(0);
        if (cost) {
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
                        (*cost)[i] *= (1._rt - 2._rt/load_balance_intervals.localPeriod(step+1));
                    }
                }
            }
        }

        // At the beginning, we have B^{n} and E^{n}.
        // Particles have p^{n} and x^{n}.
        // is_synchronized is true.
        if (is_synchronized) {
            if (do_electrostatic == ElectrostaticSolverAlgo::None) {
                // Not called at each iteration, so exchange all guard cells
                FillBoundaryE(guard_cells.ng_alloc_EB);
                FillBoundaryB(guard_cells.ng_alloc_EB);
                UpdateAuxilaryData();
                FillBoundaryAux(guard_cells.ng_UpdateAux);
            }
            // on first step, push p by -0.5*dt
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                mypc->PushP(lev, -0.5_rt*dt[lev],
                            *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                            *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
            }
            is_synchronized = false;
        } else {
            if (do_electrostatic == ElectrostaticSolverAlgo::None) {
                // Beyond one step, we have E^{n} and B^{n}.
                // Particles have p^{n-1/2} and x^{n}.

                // E and B are up-to-date inside the domain only
                FillBoundaryE(guard_cells.ng_FieldGather);
                FillBoundaryB(guard_cells.ng_FieldGather);
                // E and B: enough guard cells to update Aux or call Field Gather in fp and cp
                // Need to update Aux on lower levels, to interpolate to higher levels.
                if (fft_do_time_averaging)
                {
                    FillBoundaryE_avg(guard_cells.ng_FieldGather);
                    FillBoundaryB_avg(guard_cells.ng_FieldGather);
                }
                // TODO Remove call to FillBoundaryAux before UpdateAuxilaryData?
                if (WarpX::maxwell_solver_id != MaxwellSolverAlgo::PSATD)
                    FillBoundaryAux(guard_cells.ng_UpdateAux);
                UpdateAuxilaryData();
                FillBoundaryAux(guard_cells.ng_UpdateAux);
            }
        }

        // Run multi-physics modules:
        // ionization, Coulomb collisions, QED
        doFieldIonization();
        ExecutePythonCallback("beforecollisions");
        mypc->doCollisions( cur_time, dt[0] );
        ExecutePythonCallback("aftercollisions");
#ifdef WARPX_QED
        doQEDEvents();
        mypc->doQEDSchwinger();
#endif

        // Main PIC operation:
        // gather fields, push particles, deposit sources, update fields

        ExecutePythonCallback("particleinjection");
        // Electrostatic case: only gather fields and push particles,
        // deposition and calculation of fields done further below
        if (do_electrostatic != ElectrostaticSolverAlgo::None)
        {
            const bool skip_deposition = true;
            PushParticlesandDepose(cur_time, skip_deposition);
        }
        // Electromagnetic case: multi-J algorithm
        else if (do_multi_J)
        {
            OneStep_multiJ(cur_time);
        }
        // Electromagnetic case: no subcycling or no mesh refinement
        else if (do_subcycling == 0 || finest_level == 0)
        {
            OneStep_nosub(cur_time);
            // E: guard cells are up-to-date
            // B: guard cells are NOT up-to-date
            // F: guard cells are NOT up-to-date
        }
        // Electromagnetic case: subcycling with one level of mesh refinement
        else if (do_subcycling == 1 && finest_level == 1)
        {
            OneStep_sub1(cur_time);
        }
        else
        {
            amrex::Abort(Utils::TextMsg::Err(
                "do_subcycling = " + std::to_string(do_subcycling)
                + " is an unsupported do_subcycling type."));
        }

        // Resample particles
        // +1 is necessary here because value of step seen by user (first step is 1) is different than
        // value of step in code (first step is 0)
        mypc->doResampling(istep[0]+1);

        if (num_mirrors>0){
            applyMirrors(cur_time);
            // E : guard cells are NOT up-to-date
            // B : guard cells are NOT up-to-date
        }

        if (cur_time + dt[0] >= stop_time - 1.e-3*dt[0] || step == numsteps_max-1) {
            // At the end of last step, push p by 0.5*dt to synchronize
            FillBoundaryE(guard_cells.ng_FieldGather);
            FillBoundaryB(guard_cells.ng_FieldGather);
            if (fft_do_time_averaging)
            {
                FillBoundaryE_avg(guard_cells.ng_FieldGather);
                FillBoundaryB_avg(guard_cells.ng_FieldGather);
            }
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


        // sync up time
        for (int i = 0; i <= max_level; ++i) {
            t_new[i] = cur_time;
        }
        multi_diags->FilterComputePackFlush( step, false, true );

        bool move_j = is_synchronized;
        // If is_synchronized we need to shift j too so that next step we can evolve E by dt/2.
        // We might need to move j because we are going to make a plotfile.
        int num_moved = MoveWindow(step+1, move_j);

        mypc->ContinuousFluxInjection(cur_time, dt[0]);

        mypc->ApplyBoundaryConditions();

        // interact the particles with EB walls (if present)
#ifdef AMREX_USE_EB
        mypc->ScrapeParticles(amrex::GetVecOfConstPtrs(m_distance_to_eb));
#endif

        m_particle_boundary_buffer->gatherParticles(*mypc, amrex::GetVecOfConstPtrs(m_distance_to_eb));

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
            if (verbose) {
                amrex::Print() << Utils::TextMsg::Info("re-sorting particles");
            }
            mypc->SortParticlesByBin(sort_bin_size);
        }

        if( do_electrostatic != ElectrostaticSolverAlgo::None ) {
            ExecutePythonCallback("beforeEsolve");
            // Electrostatic solver:
            // For each species: deposit charge and add the associated space-charge
            // E and B field to the grid ; this is done at the end of the PIC
            // loop (i.e. immediately after a `Redistribute` and before particle
            // positions are next pushed) so that the particles do not deposit out of bounds
            // and so that the fields are at the correct time in the output.
            bool const reset_fields = true;
            ComputeSpaceChargeField( reset_fields );
            ExecutePythonCallback("afterEsolve");
        }

        // afterstep callback runs with the updated global time. It is included
        // in the evolve timing.
        ExecutePythonCallback("afterstep");

        /// reduced diags
        if (reduced_diags->m_plot_rd != 0)
        {
            reduced_diags->LoadBalance();
            reduced_diags->ComputeDiags(step);
            reduced_diags->WriteToFile(step);
        }
        multi_diags->FilterComputePackFlush( step );

        // inputs: unused parameters (e.g. typos) check after step 1 has finished
        if (!early_params_checked) {
            amrex::Print() << "\n"; // better: conditional \n based on return value
            amrex::ParmParse().QueryUnusedInputs();
            this->PrintGlobalWarnings("FIRST STEP"); //Print the warning list right after the first step.
            early_params_checked = true;
        }

        // create ending time stamp for calculating elapsed time each iteration
        Real evolve_time_end_step = amrex::second();
        evolve_time += evolve_time_end_step - evolve_time_beg_step;

        if (verbose) {
            amrex::Print()<< "STEP " << step+1 << " ends." << " TIME = " << cur_time
                        << " DT = " << dt[0] << "\n";
            amrex::Print()<< "Evolve time = " << evolve_time
                      << " s; This step = " << evolve_time_end_step-evolve_time_beg_step
                      << " s; Avg. per step = " << evolve_time/(step-step_begin+1) << " s\n";
        }

        if (cur_time >= stop_time - 1.e-3*dt[0]) {
            break;
        }

        // End loop on time steps
    }
    multi_diags->FilterComputePackFlushLastTimestep( istep[0] );

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
    WARPX_PROFILE("WarpX::OneStep_nosub()");

    // Push particle from x^{n} to x^{n+1}
    //               from p^{n-1/2} to p^{n+1/2}
    // Deposit current j^{n+1/2}
    // Deposit charge density rho^{n}

    ExecutePythonCallback("particlescraper");
    ExecutePythonCallback("beforedeposition");

    PushParticlesandDepose(cur_time);

    ExecutePythonCallback("afterdeposition");

    // Synchronize J and rho
    SyncCurrent();
    SyncRho();

    // At this point, J is up-to-date inside the domain, and E and B are
    // up-to-date including enough guard cells for first step of the field
    // solve.

    // For extended PML: copy J from regular grid to PML, and damp J in PML
    if (do_pml && pml_has_particles) CopyJPML();
    if (do_pml && do_pml_j_damping) DampJPML();

    ExecutePythonCallback("beforeEsolve");

    // Push E and B from {n} to {n+1}
    // (And update guard cells immediately afterwards)
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        if (use_hybrid_QED)
        {
            WarpX::Hybrid_QED_Push(dt);
            FillBoundaryE(guard_cells.ng_alloc_EB);
        }
        PushPSATD();

        if (do_pml) {
            DampPML();
        }

        if (use_hybrid_QED) {
            FillBoundaryE(guard_cells.ng_alloc_EB);
            FillBoundaryB(guard_cells.ng_alloc_EB);
            WarpX::Hybrid_QED_Push(dt);
            FillBoundaryE(guard_cells.ng_afterPushPSATD);
        }
        else {
            FillBoundaryE(guard_cells.ng_afterPushPSATD);
            FillBoundaryB(guard_cells.ng_afterPushPSATD);
            if (WarpX::do_dive_cleaning || WarpX::do_pml_dive_cleaning)
                FillBoundaryF(guard_cells.ng_alloc_F);
            if (WarpX::do_divb_cleaning || WarpX::do_pml_divb_cleaning)
                FillBoundaryG(guard_cells.ng_alloc_G);
        }

        // Synchronize E, B, F, G fields on nodal points
        NodalSync(Efield_fp, Efield_cp);
        NodalSync(Bfield_fp, Bfield_cp);
        if (WarpX::do_dive_cleaning) NodalSync(F_fp, F_cp);
        if (WarpX::do_divb_cleaning) NodalSync(G_fp, G_cp);

        if (do_pml) {
            NodalSyncPML();
        }
    } else {
        EvolveF(0.5_rt * dt[0], DtType::FirstHalf);
        EvolveG(0.5_rt * dt[0], DtType::FirstHalf);
        FillBoundaryF(guard_cells.ng_FieldSolverF);
        FillBoundaryG(guard_cells.ng_FieldSolverG);
        EvolveB(0.5_rt * dt[0], DtType::FirstHalf); // We now have B^{n+1/2}

        FillBoundaryB(guard_cells.ng_FieldSolver);

        if (WarpX::em_solver_medium == MediumForEM::Vacuum) {
            // vacuum medium
            EvolveE(dt[0]); // We now have E^{n+1}
        } else if (WarpX::em_solver_medium == MediumForEM::Macroscopic) {
            // macroscopic medium
            MacroscopicEvolveE(dt[0]); // We now have E^{n+1}
        } else {
            amrex::Abort(" Medium for EM is unknown \n");
        }

        FillBoundaryE(guard_cells.ng_FieldSolver);
        EvolveF(0.5_rt * dt[0], DtType::SecondHalf);
        EvolveG(0.5_rt * dt[0], DtType::SecondHalf);
        EvolveB(0.5_rt * dt[0], DtType::SecondHalf); // We now have B^{n+1}

        // Synchronize E and B fields on nodal points
        NodalSync(Efield_fp, Efield_cp);
        NodalSync(Bfield_fp, Bfield_cp);

        if (do_pml) {
            FillBoundaryF(guard_cells.ng_alloc_F);
            DampPML();
            NodalSyncPML();
            FillBoundaryE(guard_cells.ng_MovingWindow);
            FillBoundaryF(guard_cells.ng_MovingWindow);
            FillBoundaryB(guard_cells.ng_MovingWindow);
        }
        // E and B are up-to-date in the domain, but all guard cells are
        // outdated.
        if (safe_guard_cells)
            FillBoundaryB(guard_cells.ng_alloc_EB);
    } // !PSATD

    ExecutePythonCallback("afterEsolve");
}

void
WarpX::OneStep_multiJ (const amrex::Real cur_time)
{
#ifdef WARPX_USE_PSATD
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD)
    {
        // Push particle from x^{n} to x^{n+1}
        //               from p^{n-1/2} to p^{n+1/2}
        const bool skip_deposition = true;
        PushParticlesandDepose(cur_time, skip_deposition);

        // Initialize multi-J loop:

        // 1) Prepare E,B,F,G fields in spectral space
        PSATDForwardTransformEB();
        if (WarpX::do_dive_cleaning) PSATDForwardTransformF();
        if (WarpX::do_divb_cleaning) PSATDForwardTransformG();

        // 2) Set the averaged fields to zero
        if (WarpX::fft_do_time_averaging) PSATDEraseAverageFields();

        // 3) Deposit rho (in rho_new, since it will be moved during the loop)
        if (WarpX::update_with_rho)
        {
            // Deposit rho at relative time -dt
            // (dt[0] denotes the time step on mesh refinement level 0)
            mypc->DepositCharge(rho_fp, -dt[0]);
            // Filter, exchange boundary, and interpolate across levels
            SyncRho();
            // Forward FFT of rho_new
            PSATDForwardTransformRho(0, 1);
        }

        // 4) Deposit J at relative time -dt with time step dt
        //    (dt[0] denotes the time step on mesh refinement level 0)
        auto& current = (WarpX::do_current_centering) ? current_fp_nodal : current_fp;
        mypc->DepositCurrent(current, dt[0], -dt[0]);
        // Filter, exchange boundary, and interpolate across levels
        SyncCurrent();
        // Forward FFT of J
        PSATDForwardTransformJ(current_fp, current_cp);

        // Number of depositions for multi-J scheme
        const int n_depose = WarpX::do_multi_J_n_depositions;
        // Time sub-step for each multi-J deposition
        const amrex::Real sub_dt = dt[0] / static_cast<amrex::Real>(n_depose);
        // Whether to perform multi-J depositions on a time interval that spans
        // one or two full time steps (from n*dt to (n+1)*dt, or from n*dt to (n+2)*dt)
        const int n_loop = (WarpX::fft_do_time_averaging) ? 2*n_depose : n_depose;

        // Loop over multi-J depositions
        for (int i_depose = 0; i_depose < n_loop; i_depose++)
        {
            // Move J deposited previously, from new to old
            PSATDMoveJNewToJOld();

            const amrex::Real t_depose = (i_depose-n_depose+1)*sub_dt;

            // Deposit new J at relative time t_depose with time step dt
            // (dt[0] denotes the time step on mesh refinement level 0)
            mypc->DepositCurrent(current, dt[0], t_depose);
            // Filter, exchange boundary, and interpolate across levels
            SyncCurrent();
            // Forward FFT of J
            PSATDForwardTransformJ(current_fp, current_cp);

            // Deposit new rho
            if (WarpX::update_with_rho)
            {
                // Move rho deposited previously, from new to old
                PSATDMoveRhoNewToRhoOld();

                // Deposit rho at relative time t_depose
                mypc->DepositCharge(rho_fp, t_depose);
                // Filter, exchange boundary, and interpolate across levels
                SyncRho();
                // Forward FFT of rho_new
                PSATDForwardTransformRho(0, 1);
            }

            // Advance E,B,F,G fields in time and update the average fields
            PSATDPushSpectralFields();

            // Transform non-average fields E,B,F,G after n_depose pushes
            // (the relative time reached here coincides with an integer full time step)
            if (i_depose == n_depose-1)
            {
                PSATDBackwardTransformEB();
                if (WarpX::do_dive_cleaning) PSATDBackwardTransformF();
                if (WarpX::do_divb_cleaning) PSATDBackwardTransformG();
            }
        }

        // Transform fields back to real space
        if (WarpX::fft_do_time_averaging)
        {
            // We summed the integral of the field over 2*dt
            PSATDScaleAverageFields(1._rt / (2._rt*dt[0]));
            PSATDBackwardTransformEBavg();
        }

        // Evolve fields in PML
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (do_pml && pml[lev]->ok())
            {
                pml[lev]->PushPSATD(lev);
            }
            ApplyEfieldBoundary(lev, PatchType::fine);
            if (lev > 0) ApplyEfieldBoundary(lev, PatchType::coarse);
            ApplyBfieldBoundary(lev, PatchType::fine, DtType::FirstHalf);
            if (lev > 0) ApplyBfieldBoundary(lev, PatchType::coarse, DtType::FirstHalf);
        }

        // Damp fields in PML before exchanging guard cells
        if (do_pml)
        {
            DampPML();
        }

        // Exchange guard cells
        FillBoundaryE(guard_cells.ng_alloc_EB);
        FillBoundaryB(guard_cells.ng_alloc_EB);
        if (WarpX::do_dive_cleaning || WarpX::do_pml_dive_cleaning) FillBoundaryF(guard_cells.ng_alloc_F);
        if (WarpX::do_divb_cleaning || WarpX::do_pml_divb_cleaning) FillBoundaryG(guard_cells.ng_alloc_G);

        // Synchronize E, B, F, G fields on nodal points
        NodalSync(Efield_fp, Efield_cp);
        NodalSync(Bfield_fp, Bfield_cp);
        if (WarpX::do_dive_cleaning) NodalSync(F_fp, F_cp);
        if (WarpX::do_divb_cleaning) NodalSync(G_fp, G_cp);

        // Synchronize fields on nodal points in PML
        if (do_pml)
        {
            NodalSyncPML();
        }
    }
    else
    {
        amrex::Abort("multi-J algorithm not implemented for FDTD");
    }
#else
    amrex::ignore_unused(cur_time);
    amrex::Abort("multi-J algorithm not implemented for FDTD");
#endif // WARPX_USE_PSATD
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

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(finest_level == 1, "Must have exactly two levels");
    const int fine_lev = 1;
    const int coarse_lev = 0;

    // i) Push particles and fields on the fine patch (first fine step)
    PushParticlesandDepose(fine_lev, curtime, DtType::FirstHalf);
    RestrictCurrentFromFineToCoarsePatch(fine_lev);
    RestrictRhoFromFineToCoarsePatch(fine_lev);
    ApplyFilterandSumBoundaryJ(current_fp, current_cp, fine_lev, PatchType::fine);
    NodalSyncJ(current_fp, current_cp, fine_lev, PatchType::fine);
    ApplyFilterandSumBoundaryRho(fine_lev, PatchType::fine, 0, 2*ncomps);
    NodalSyncRho(fine_lev, PatchType::fine, 0, 2);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::FirstHalf);
    EvolveF(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);
    FillBoundaryF(fine_lev, PatchType::fine, guard_cells.ng_alloc_F);

    EvolveE(fine_lev, PatchType::fine, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::fine, guard_cells.ng_FieldGather);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::SecondHalf);
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
    AddCurrentFromFineLevelandSumBoundary(current_fp, current_cp, coarse_lev);
    AddRhoFromFineLevelandSumBoundary(coarse_lev, 0, ncomps);

    EvolveB(fine_lev, PatchType::coarse, dt[fine_lev], DtType::FirstHalf);
    EvolveF(fine_lev, PatchType::coarse, dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::coarse, guard_cells.ng_FieldGather);
    FillBoundaryF(fine_lev, PatchType::coarse, guard_cells.ng_FieldSolverF);

    EvolveE(fine_lev, PatchType::coarse, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::coarse, guard_cells.ng_FieldGather);

    EvolveB(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev], DtType::FirstHalf);
    EvolveF(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev], DtType::FirstHalf);
    FillBoundaryB(coarse_lev, PatchType::fine, guard_cells.ng_FieldGather);
    FillBoundaryF(coarse_lev, PatchType::fine, guard_cells.ng_FieldSolverF);

    EvolveE(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev]);
    FillBoundaryE(coarse_lev, PatchType::fine, guard_cells.ng_FieldGather);

    // TODO Remove call to FillBoundaryAux before UpdateAuxilaryData?
    FillBoundaryAux(guard_cells.ng_UpdateAux);
    // iii) Get auxiliary fields on the fine grid, at dt[fine_lev]
    UpdateAuxilaryData();
    FillBoundaryAux(guard_cells.ng_UpdateAux);

    // iv) Push particles and fields on the fine patch (second fine step)
    PushParticlesandDepose(fine_lev, curtime+dt[fine_lev], DtType::SecondHalf);
    RestrictCurrentFromFineToCoarsePatch(fine_lev);
    RestrictRhoFromFineToCoarsePatch(fine_lev);
    ApplyFilterandSumBoundaryJ(current_fp, current_cp, fine_lev, PatchType::fine);
    NodalSyncJ(current_fp, current_cp, fine_lev, PatchType::fine);
    ApplyFilterandSumBoundaryRho(fine_lev, PatchType::fine, 0, ncomps);
    NodalSyncRho(fine_lev, PatchType::fine, 0, 2);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::FirstHalf);
    EvolveF(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);
    FillBoundaryF(fine_lev, PatchType::fine, guard_cells.ng_FieldSolverF);

    EvolveE(fine_lev, PatchType::fine, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::fine, guard_cells.ng_FieldSolver);

    EvolveB(fine_lev, PatchType::fine, 0.5_rt*dt[fine_lev], DtType::SecondHalf);
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
    AddCurrentFromFineLevelandSumBoundary(current_fp, current_cp, coarse_lev);
    AddRhoFromFineLevelandSumBoundary(coarse_lev, ncomps, ncomps);

    EvolveE(fine_lev, PatchType::coarse, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::coarse, guard_cells.ng_FieldSolver);

    EvolveB(fine_lev, PatchType::coarse, dt[fine_lev], DtType::SecondHalf);
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

    EvolveB(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev], DtType::SecondHalf);
    EvolveF(coarse_lev, PatchType::fine, 0.5_rt*dt[coarse_lev], DtType::SecondHalf);

    if (do_pml) {
        if (moving_window_active(istep[0]+1)){
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

    // Synchronize all nodal points at the end of the timestep
    NodalSync(Efield_fp, Efield_cp);
    NodalSync(Bfield_fp, Bfield_cp);
    if (WarpX::do_dive_cleaning) NodalSync(F_fp, F_cp);
    if (do_pml) NodalSyncPML();
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
WarpX::PushParticlesandDepose (amrex::Real cur_time, bool skip_deposition)
{
    // Evolve particles to p^{n+1/2} and x^{n+1}
    // Depose current, j^{n+1/2}
    for (int lev = 0; lev <= finest_level; ++lev) {
        PushParticlesandDepose(lev, cur_time, DtType::Full, skip_deposition);
    }
}

void
WarpX::PushParticlesandDepose (int lev, amrex::Real cur_time, DtType a_dt_type, bool skip_deposition)
{
    amrex::MultiFab* current_x = nullptr;
    amrex::MultiFab* current_y = nullptr;
    amrex::MultiFab* current_z = nullptr;

    if (WarpX::do_current_centering)
    {
        current_x = current_fp_nodal[lev][0].get();
        current_y = current_fp_nodal[lev][1].get();
        current_z = current_fp_nodal[lev][2].get();
    }
    else if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Vay)
    {
        current_x = current_fp_vay[lev][0].get();
        current_y = current_fp_vay[lev][1].get();
        current_z = current_fp_vay[lev][2].get();
    }
    else
    {
        current_x = current_fp[lev][0].get();
        current_y = current_fp[lev][1].get();
        current_z = current_fp[lev][2].get();
    }

    mypc->Evolve(lev,
                 *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                 *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2],
                 *current_x, *current_y, *current_z,
                 current_buf[lev][0].get(), current_buf[lev][1].get(), current_buf[lev][2].get(),
                 rho_fp[lev].get(), charge_buf[lev].get(),
                 Efield_cax[lev][0].get(), Efield_cax[lev][1].get(), Efield_cax[lev][2].get(),
                 Bfield_cax[lev][0].get(), Bfield_cax[lev][1].get(), Bfield_cax[lev][2].get(),
                 cur_time, dt[lev], a_dt_type, skip_deposition);
#ifdef WARPX_DIM_RZ
    if (! skip_deposition) {
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
    }
#endif
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
