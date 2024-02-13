/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "Diagnostics/MultiDiagnostics.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "Evolve/WarpXDtType.H"
#include "Evolve/WarpXPushType.H"
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
#include "Python/callbacks.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

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

void
WarpX::EvolveImplicitEMInit (int lev)
{

    if (lev == 0) {
        // Add space to save the positions and velocities at the start of the time steps
        for (auto const& pc : *mypc) {
#if (AMREX_SPACEDIM >= 2)
            pc->AddRealComp("x_n");
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
            pc->AddRealComp("y_n");
#endif
            pc->AddRealComp("z_n");
            pc->AddRealComp("ux_n");
            pc->AddRealComp("uy_n");
            pc->AddRealComp("uz_n");
        }
    }

    // Initialize MultiFabs to hold the E and B fields at the start of the time steps
    // Only one refinement level is supported
    const int nlevs_max = maxLevel() + 1;
    Efield_rhs.resize(nlevs_max);
    Efield_n.resize(nlevs_max);
    Efield_save.resize(nlevs_max);
    Bfield_rhs.resize(nlevs_max);
    if (evolve_scheme == EvolveScheme::ThetaImplicit) {
        Bfield_n.resize(nlevs_max);
        Bfield_save.resize(nlevs_max);
    }

    // The Efield_n and Bfield_n will hold the fields at the start of the time step.
    // This is needed since in each iteration the fields are advanced from the values
    // at the start of the step.
    // The Efield_save and Bfield_save will hold the fields from the previous iteration,
    // to check the change in the fields after the iterations to check for convergence.
    // The Efiel_fp and Bfield_fp will hole the n+theta during the iterations and then
    // advance to the n+1 time level after the iterations complete.
    AllocInitMultiFabFromModel(Efield_n[lev][0], *Efield_fp[0][0], lev, "Efield_n[0]");
    AllocInitMultiFabFromModel(Efield_n[lev][1], *Efield_fp[0][1], lev, "Efield_n[1]");
    AllocInitMultiFabFromModel(Efield_n[lev][2], *Efield_fp[0][2], lev, "Efield_n[2]");
    AllocInitMultiFabFromModel(Efield_rhs[lev][0], *Efield_fp[0][0], lev, "Efield_rhs[0]");
    AllocInitMultiFabFromModel(Efield_rhs[lev][1], *Efield_fp[0][1], lev, "Efield_rhs[1]");
    AllocInitMultiFabFromModel(Efield_rhs[lev][2], *Efield_fp[0][2], lev, "Efield_rhs[2]");
    AllocInitMultiFabFromModel(Efield_save[lev][0], *Efield_fp[0][0], lev, "Efield_save[0]");
    AllocInitMultiFabFromModel(Efield_save[lev][1], *Efield_fp[0][1], lev, "Efield_save[1]");
    AllocInitMultiFabFromModel(Efield_save[lev][2], *Efield_fp[0][2], lev, "Efield_save[2]");

    if (evolve_scheme == EvolveScheme::ThetaImplicit) {
        AllocInitMultiFabFromModel(Bfield_n[lev][0], *Bfield_fp[0][0], lev, "Bfield_n[0]");
        AllocInitMultiFabFromModel(Bfield_n[lev][1], *Bfield_fp[0][1], lev, "Bfield_n[1]");
        AllocInitMultiFabFromModel(Bfield_n[lev][2], *Bfield_fp[0][2], lev, "Bfield_n[2]");
        AllocInitMultiFabFromModel(Bfield_save[lev][0], *Bfield_fp[0][0], lev, "Bfield_save[0]");
        AllocInitMultiFabFromModel(Bfield_save[lev][1], *Bfield_fp[0][1], lev, "Bfield_save[1]");
        AllocInitMultiFabFromModel(Bfield_save[lev][2], *Bfield_fp[0][2], lev, "Bfield_save[2]");
    }
    AllocInitMultiFabFromModel(Bfield_rhs[lev][0], *Bfield_fp[0][0], lev, "Bfield_rhs[0]");
    AllocInitMultiFabFromModel(Bfield_rhs[lev][1], *Bfield_fp[0][1], lev, "Bfield_rhs[1]");
    AllocInitMultiFabFromModel(Bfield_rhs[lev][2], *Bfield_fp[0][2], lev, "Bfield_rhs[2]");

}

void
WarpX::OneStep_ImplicitEM(amrex::Real cur_time)
{
    using namespace amrex::literals;

    // We have E^{n}.
    // Particles have p^{n} and x^{n}.
    // With full implicit, B^{n}
    // With semi-implicit, B^{n-1/2}

    // Save the values at the start of the time step,
    // copying particle data to x_n etc.
    for (auto const& pc : *mypc) {
        SaveParticlesAtImplicitStepStart (*pc, 0);
    }

    // Save the fields at the start of the step
    amrex::MultiFab::Copy(*Efield_n[0][0], *Efield_fp[0][0], 0, 0, ncomps, Efield_fp[0][0]->nGrowVect());
    amrex::MultiFab::Copy(*Efield_n[0][1], *Efield_fp[0][1], 0, 0, ncomps, Efield_fp[0][1]->nGrowVect());
    amrex::MultiFab::Copy(*Efield_n[0][2], *Efield_fp[0][2], 0, 0, ncomps, Efield_fp[0][2]->nGrowVect());

    if (evolve_scheme == EvolveScheme::ThetaImplicit) {
        amrex::MultiFab::Copy(*Bfield_n[0][0], *Bfield_fp[0][0], 0, 0, ncomps, Bfield_fp[0][0]->nGrowVect());
        amrex::MultiFab::Copy(*Bfield_n[0][1], *Bfield_fp[0][1], 0, 0, ncomps, Bfield_fp[0][1]->nGrowVect());
        amrex::MultiFab::Copy(*Bfield_n[0][2], *Bfield_fp[0][2], 0, 0, ncomps, Bfield_fp[0][2]->nGrowVect());
    } else if (evolve_scheme == EvolveScheme::SemiImplicit) {
        // Compute Bfield at time n+1/2
        ComputeRHSB(dt[0]);
        Bfield_fp[0][0]->plus(*Bfield_rhs[0][0], 0, ncomps, 0);
        Bfield_fp[0][1]->plus(*Bfield_rhs[0][1], 0, ncomps, 0);
        Bfield_fp[0][2]->plus(*Bfield_rhs[0][2], 0, ncomps, 0);

        // WarpX::sync_nodal_points is used to avoid instability
        FillBoundaryB(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
        ApplyBfieldBoundary(0, PatchType::fine, DtType::Full);
    }

    // Start the iterations
    amrex::Real deltaE = 1._rt;
    amrex::Real deltaB = 1._rt;
    int iteration_count = 0;
    while (iteration_count < max_picard_iterations &&
           (deltaE > picard_iteration_tolerance || deltaB > picard_iteration_tolerance)) {
        iteration_count++;

        // Advance the particle positions by 1/2 dt,
        // particle velocities by dt, then take average of old and new v,
        // deposit currents, giving J at n+1/2 used in ComputeRHSE below
        PreRHSOpFromNonlinearIter( cur_time, dt[0], iteration_count );

        if (picard_iteration_tolerance > 0. || iteration_count == max_picard_iterations) {
            // Save the E at n+1/2 from the previous iteration so that the change
            // in this iteration can be calculated
            amrex::MultiFab::Copy(*Efield_save[0][0], *Efield_fp[0][0], 0, 0, ncomps, 0);
            amrex::MultiFab::Copy(*Efield_save[0][1], *Efield_fp[0][1], 0, 0, ncomps, 0);
            amrex::MultiFab::Copy(*Efield_save[0][2], *Efield_fp[0][2], 0, 0, ncomps, 0);
        }

        // Compute Efield at time n+1/2
        amrex::MultiFab::Copy(*Efield_fp[0][0], *Efield_n[0][0], 0, 0, ncomps, Efield_n[0][0]->nGrowVect());
        amrex::MultiFab::Copy(*Efield_fp[0][1], *Efield_n[0][1], 0, 0, ncomps, Efield_n[0][1]->nGrowVect());
        amrex::MultiFab::Copy(*Efield_fp[0][2], *Efield_n[0][2], 0, 0, ncomps, Efield_n[0][2]->nGrowVect());

        ComputeRHSE(0.5_rt*dt[0]);
        Efield_fp[0][0]->plus(*Efield_rhs[0][0], 0, ncomps, 0);
        Efield_fp[0][1]->plus(*Efield_rhs[0][1], 0, ncomps, 0);
        Efield_fp[0][2]->plus(*Efield_rhs[0][2], 0, ncomps, 0);

        // WarpX::sync_nodal_points is used to avoid instability
        FillBoundaryE(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
        ApplyEfieldBoundary(0, PatchType::fine);

        if (evolve_scheme == EvolveScheme::ThetaImplicit) {
            if (picard_iteration_tolerance > 0. || iteration_count == max_picard_iterations) {
                // Save the B at n+1/2 from the previous iteration so that the change
                // in this iteration can be calculated
                amrex::MultiFab::Copy(*Bfield_save[0][0], *Bfield_fp[0][0], 0, 0, ncomps, 0);
                amrex::MultiFab::Copy(*Bfield_save[0][1], *Bfield_fp[0][1], 0, 0, ncomps, 0);
                amrex::MultiFab::Copy(*Bfield_save[0][2], *Bfield_fp[0][2], 0, 0, ncomps, 0);
            }

            // Compute Bfield at time n+1/2
            amrex::MultiFab::Copy(*Bfield_fp[0][0], *Bfield_n[0][0], 0, 0, ncomps, Bfield_n[0][0]->nGrowVect());
            amrex::MultiFab::Copy(*Bfield_fp[0][1], *Bfield_n[0][1], 0, 0, ncomps, Bfield_n[0][1]->nGrowVect());
            amrex::MultiFab::Copy(*Bfield_fp[0][2], *Bfield_n[0][2], 0, 0, ncomps, Bfield_n[0][2]->nGrowVect());

            ComputeRHSB(0.5_rt*dt[0]);
            Bfield_fp[0][0]->plus(*Bfield_rhs[0][0], 0, ncomps, 0);
            Bfield_fp[0][1]->plus(*Bfield_rhs[0][1], 0, ncomps, 0);
            Bfield_fp[0][2]->plus(*Bfield_rhs[0][2], 0, ncomps, 0);

            // WarpX::sync_nodal_points is used to avoid instability
            FillBoundaryB(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
            ApplyBfieldBoundary(0, PatchType::fine, DtType::Full);
        }

        // The B field update needs
        if (num_mirrors>0){
            applyMirrors(cur_time);
            // E : guard cells are NOT up-to-date from the mirrors
            // B : guard cells are NOT up-to-date from the mirrors
        }

        if (picard_iteration_tolerance > 0. || iteration_count == max_picard_iterations) {
            // Calculate the change in E and B from this iteration
            // deltaE = abs(Enew - Eold)/max(abs(Enew))
            Efield_save[0][0]->minus(*Efield_fp[0][0], 0, ncomps, 0);
            Efield_save[0][1]->minus(*Efield_fp[0][1], 0, ncomps, 0);
            Efield_save[0][2]->minus(*Efield_fp[0][2], 0, ncomps, 0);
            amrex::Real maxE0 = std::max(1._rt, Efield_fp[0][0]->norm0(0, 0));
            amrex::Real maxE1 = std::max(1._rt, Efield_fp[0][1]->norm0(0, 0));
            amrex::Real maxE2 = std::max(1._rt, Efield_fp[0][2]->norm0(0, 0));
            amrex::Real deltaE0 = Efield_save[0][0]->norm0(0, 0)/maxE0;
            amrex::Real deltaE1 = Efield_save[0][1]->norm0(0, 0)/maxE1;
            amrex::Real deltaE2 = Efield_save[0][2]->norm0(0, 0)/maxE2;
            deltaE = std::max(std::max(deltaE0, deltaE1), deltaE2);
            if (evolve_scheme == EvolveScheme::ThetaImplicit) {
                Bfield_save[0][0]->minus(*Bfield_fp[0][0], 0, ncomps, 0);
                Bfield_save[0][1]->minus(*Bfield_fp[0][1], 0, ncomps, 0);
                Bfield_save[0][2]->minus(*Bfield_fp[0][2], 0, ncomps, 0);
                amrex::Real maxB0 = std::max(1._rt, Bfield_fp[0][0]->norm0(0, 0));
                amrex::Real maxB1 = std::max(1._rt, Bfield_fp[0][1]->norm0(0, 0));
                amrex::Real maxB2 = std::max(1._rt, Bfield_fp[0][2]->norm0(0, 0));
                amrex::Real deltaB0 = Bfield_save[0][0]->norm0(0, 0)/maxB0;
                amrex::Real deltaB1 = Bfield_save[0][1]->norm0(0, 0)/maxB1;
                amrex::Real deltaB2 = Bfield_save[0][2]->norm0(0, 0)/maxB2;
                deltaB = std::max(std::max(deltaB0, deltaB1), deltaB2);
            } else {
                deltaB = 0.;
            }
            amrex::Print() << "Max delta " << iteration_count << " " << deltaE << " " << deltaB << "\n";
        }

        // Now, the particle positions and velocities and the Efield_fp and Bfield_fp hold
        // the new values at n+1/2
    }

    amrex::Print() << "Picard iterations = " << iteration_count << ", Eerror = " << deltaE << ", Berror = " << deltaB << "\n";
    if (picard_iteration_tolerance > 0. && iteration_count == max_picard_iterations) {
       std::stringstream convergenceMsg;
       convergenceMsg << "The Picard implicit solver failed to converge after " << iteration_count << " iterations, with Eerror = " << deltaE << ", Berror = " << deltaB << " with a tolerance of " << picard_iteration_tolerance;
       if (require_picard_convergence) {
           WARPX_ABORT_WITH_MESSAGE(convergenceMsg.str());
       } else {
           ablastr::warn_manager::WMRecordWarning("PicardSolver", convergenceMsg.str());
       }
    }

    // Advance particles to step n+1
    for (auto const& pc : *mypc) {
        FinishImplicitParticleUpdate(*pc, 0);
    }

    // Advance fields to step n+1
    // WarpX::sync_nodal_points is used to avoid instability
    FinishImplicitFieldUpdate(Efield_fp, Efield_n);
    FillBoundaryE(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
    if (evolve_scheme == EvolveScheme::ThetaImplicit) {
        FinishImplicitFieldUpdate(Bfield_fp, Bfield_n);
        FillBoundaryB(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
    }

}
