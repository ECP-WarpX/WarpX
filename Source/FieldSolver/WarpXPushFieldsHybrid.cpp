/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Evolve/WarpXDtType.H"
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"

using namespace amrex;

void
WarpX::HybridEvolveFields ()
{
    // HybridSolveE(DtType::FirstHalf);
    // FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // EvolveB(dt[0], DtType::FirstHalf);
    // FillBoundaryB(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // return;

    // During the particle push and deposition (which already happened) the
    // charge density and current density was updated. So at this time we
    // have rho^{n} in the 0'th index and rho{n+1} in the 1'st index of `rho_fp`,
    // J^{n+1/2} in `current_fp` and J^{n-1/2} in `current_fp_old`.

    // Firstly, E^{n} is recalculated with the accurate V^{n} since at the end
    // of the last step we had to "guess" V^{n}.
    // HybridSolveE(DtType::Full);
    // FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // TODO: insert Runge-Kutta integration logic to supercycle B update instead
    // of the single update step used here - can test with small timestep using
    // this simpler implementation

    // Push the B field from t=n to t=n+1/2 using the current and density
    // at t=n, but updating the E field along with B.
    for (int sub_step = 0; sub_step < 5; sub_step++)
    {
        HybridSolveE(DtType::Full);
        FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

        EvolveB(0.1_rt * dt[0], DtType::FirstHalf);
        FillBoundaryB(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);
    }

    // Now push the B field to t=n+1 using the n+1/2 quantities
    for (int sub_step = 0; sub_step < 5; sub_step++)
    {
        HybridSolveE(DtType::FirstHalf);
        FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

        EvolveB(0.1_rt * dt[0], DtType::SecondHalf);
        FillBoundaryB(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);
    }

    // First the B field is pushed forward half a timestep
    // to get B^{n+1/2}
    // EvolveB(0.5_rt * dt[0], DtType::FirstHalf);
    // FillBoundaryB(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // // Now the E field is updated using the electron momentum equation
    // HybridSolveE(DtType::FirstHalf);
    // FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // // The B field is pushed forward another half a timestep
    // // to get B^{n+1}
    // EvolveB(0.5_rt * dt[0], DtType::SecondHalf);
    // FillBoundaryB(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // Update the E field to E^{n+1}
    HybridSolveE(DtType::SecondHalf);
    FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // Finally, the "new" current density values are copied to the "old"
    // location
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::Copy(*current_fp_old[lev][idim], *current_fp[lev][idim],
                           0, 0, 1, current_fp_old[lev][idim]->nGrowVect());
        }
    }
}

void
WarpX::HybridSolveE (DtType a_dt_type)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        HybridSolveE(lev, a_dt_type);
    }
}

void
WarpX::HybridSolveE (int lev, DtType a_dt_type)
{
    WARPX_PROFILE("WarpX::HybridSolveE()");
    HybridSolveE(lev, PatchType::fine, a_dt_type);
    if (lev > 0)
    {
        amrex::Abort(Utils::TextMsg::Err(
        "HybridSolveE: Only one level implemented for hybrid solver."));
        HybridSolveE(lev, PatchType::coarse, a_dt_type);
    }
}

void
WarpX::HybridSolveE (int lev, PatchType patch_type, DtType a_dt_type)
{
    // Solve E field in regular cells
    if (patch_type == PatchType::fine) {
        m_fdtd_solver_fp[lev]->HybridSolveE(
            Efield_fp[lev], Bfield_fp[lev],
            current_fp[lev], current_fp_old[lev], rho_fp[lev],
            m_edge_lengths[lev], lev, m_hybrid_model, a_dt_type
        );
    } else {
        amrex::Abort(Utils::TextMsg::Err(
        "HybridSolveE: Only one level implemented for hybrid solver."));
    }

    // Evolve E field in PML cells
    // if (do_pml && pml[lev]->ok()) {
    //     if (patch_type == PatchType::fine) {
    //         m_fdtd_solver_fp[lev]->EvolveEPML(
    //             pml[lev]->GetE_fp(), pml[lev]->GetB_fp(),
    //             pml[lev]->Getj_fp(), pml[lev]->Get_edge_lengths(),
    //             pml[lev]->GetF_fp(),
    //             pml[lev]->GetMultiSigmaBox_fp(),
    //             a_dt, pml_has_particles );
    //     } else {
    //         m_fdtd_solver_cp[lev]->EvolveEPML(
    //             pml[lev]->GetE_cp(), pml[lev]->GetB_cp(),
    //             pml[lev]->Getj_cp(), pml[lev]->Get_edge_lengths(),
    //             pml[lev]->GetF_cp(),
    //             pml[lev]->GetMultiSigmaBox_cp(),
    //             a_dt, pml_has_particles );
    //     }
    // }

    ApplyEfieldBoundary(lev, patch_type);

    // ECTRhofield must be recomputed at the very end of the Efield update to ensure
    // that ECTRhofield is consistent with Efield
// #ifdef AMREX_USE_EB
//     if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT) {
//         if (patch_type == PatchType::fine) {
//             m_fdtd_solver_fp[lev]->EvolveECTRho(Efield_fp[lev], m_edge_lengths[lev],
//                                                 m_face_areas[lev], ECTRhofield[lev], lev);
//         } else {
//             m_fdtd_solver_cp[lev]->EvolveECTRho(Efield_cp[lev], m_edge_lengths[lev],
//                                                 m_face_areas[lev], ECTRhofield[lev], lev);
//         }
//     }
// #endif
}