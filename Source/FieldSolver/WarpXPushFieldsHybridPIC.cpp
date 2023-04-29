/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Evolve/WarpXDtType.H"
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

using namespace amrex;

void WarpX::HybridPICEvolveFields ()
{
    // The below deposition is hard coded for a single level simulation
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        finest_level == 0,
        "Ohm's law E-solve only works with a single level.");

    // The particles have now been pushed to their t_{n+1} positions.
    // Perform charge deposition in component 0 of rho_fp at t_{n+1}.
    mypc->DepositCharge(rho_fp, 0._rt);
    // Perform current deposition at t_{n+1/2}.
    mypc->DepositCurrent(current_fp, dt[0], -0.5_rt * dt[0]);

    // Synchronize J and rho:
    // filter (if used), exchange guard cells, interpolate across MR levels
    // and apply boundary conditions
    SyncCurrentAndRho();

    // Get requested number of substeps to use
    int sub_steps = m_hybrid_pic_model->m_substeps / 2;

    // During the above deposition the charge and current density were updated
    // so that, at this time, we have rho^{n} in rho_fp_temp, rho{n+1} in the
    // 0'th index of `rho_fp`, J_i^{n-1/2} in `current_fp_temp` and J_i^{n+1/2}
    // in `current_fp`.

    // TODO: To speed up the algorithm insert Runge-Kutta integration logic
    // for B update instead of the substep update used here - can test with
    // small timestep using this simpler implementation

    // Note: E^{n} is recalculated with the accurate J_i^{n} since at the end
    // of the last step we had to "guess" it. It also needs to be
    // recalculated to include the resistivity before evolving B.

    // J_i^{n} is calculated as the average of J_i^{n-1/2} and J_i^{n+1/2}.
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim < 3; ++idim) {
            // Perform a linear combination of values in the 0'th index (1 comp)
            // of J_i^{n-1/2} and J_i^{n+1/2} (with 0.5 prefactors), writing
            // the result into the 0'th index of `current_fp_temp[lev][idim]`
            MultiFab::LinComb(
                *current_fp_temp[lev][idim],
                0.5_rt, *current_fp_temp[lev][idim], 0,
                0.5_rt, *current_fp[lev][idim], 0,
                0, 1, current_fp_temp[lev][idim]->nGrowVect()
            );
        }
    }

    // Calculate the electron pressure at t=n using rho^n
    CalculateElectronPressure(DtType::FirstHalf);

    // Push the B field from t=n to t=n+1/2 using the current and density
    // at t=n, while updating the E field along with B using the electron
    // momentum equation
    for (int sub_step = 0; sub_step < sub_steps; sub_step++)
    {
        CalculateCurrentAmpere();
        HybridPICSolveE(DtType::FirstHalf);
        EvolveB(0.5 / sub_steps * dt[0], DtType::FirstHalf);
        FillBoundaryB(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);
    }

    // Average rho^{n} and rho^{n+1} to get rho^{n+1/2} in rho_fp_temp
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Perform a linear combination of values in the 0'th index (1 comp)
        // of rho^{n} and rho^{n+1} (with 0.5 prefactors), writing
        // the result into the 0'th index of `rho_fp_temp[lev]`
        MultiFab::LinComb(
            *rho_fp_temp[lev], 0.5_rt, *rho_fp_temp[lev], 0,
            0.5_rt, *rho_fp[lev], 0, 0, 1, rho_fp_temp[lev]->nGrowVect()
        );
    }

    // Calculate the electron pressure at t=n+1/2
    CalculateElectronPressure(DtType::SecondHalf);

    // Now push the B field from t=n+1/2 to t=n+1 using the n+1/2 quantities
    for (int sub_step = 0; sub_step < sub_steps; sub_step++)
    {
        CalculateCurrentAmpere();
        HybridPICSolveE(DtType::SecondHalf);
        EvolveB(0.5 / sub_steps * dt[0], DtType::SecondHalf);
        FillBoundaryB(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);
    }

    // Extrapolate the ion current density to t=n+1 using
    // J_i^{n+1} = 1/2 * J_i^{n-1/2} + 3/2 * J_i^{n+1/2}, and recalling that
    // now current_fp_temp = J_i^{n} = 1/2 * (J_i^{n-1/2} + J_i^{n+1/2})
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim < 3; ++idim) {
            // Perform a linear combination of values in the 0'th index (1 comp)
            // of J_i^{n-1/2} and J_i^{n+1/2} (with -1.0 and 2.0 prefactors),
            // writing the result into the 0'th index of `current_fp_temp[lev][idim]`
            MultiFab::LinComb(
                *current_fp_temp[lev][idim],
                -1._rt, *current_fp_temp[lev][idim], 0,
                2._rt, *current_fp[lev][idim], 0,
                0, 1, current_fp_temp[lev][idim]->nGrowVect()
            );
        }
    }

    // Calculate the electron pressure at t=n+1
    CalculateElectronPressure(DtType::Full);

    // Update the E field to t=n+1 using the extrapolated J_i^n+1 value
    CalculateCurrentAmpere();
    HybridPICSolveE(DtType::Full);

    // Copy the rho^{n+1} values to rho_fp_temp and the J_i^{n+1/2} values to
    // current_fp_temp since at the next step those values will be needed as
    // rho^{n} and J_i^{n-1/2}.
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // copy 1 component value starting at index 0 to index 0
        MultiFab::Copy(*rho_fp_temp[lev], *rho_fp[lev],
                        0, 0, 1, rho_fp_temp[lev]->nGrowVect());
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::Copy(*current_fp_temp[lev][idim], *current_fp[lev][idim],
                           0, 0, 1, current_fp_temp[lev][idim]->nGrowVect());
        }
    }
}

void WarpX::CalculateCurrentAmpere ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        CalculateCurrentAmpere(lev);
    }
}

void WarpX::CalculateCurrentAmpere (const int lev)
{
    WARPX_PROFILE("WarpX::CalculateCurrentAmpere()");

    m_fdtd_solver_fp[lev]->CalculateCurrentAmpere(
        current_fp_ampere[lev], Bfield_fp[lev],
        m_edge_lengths[lev], lev
    );

    // we shouldn't apply the boundary condition to J since J = J_i - J_e but
    // the boundary correction was already applied to J_i and the B-field
    // boundary ensures that J itself complies with the boundary conditions, right?
    // ApplyJfieldBoundary(lev, Jfield[0].get(), Jfield[1].get(), Jfield[2].get());
    for (int i=0; i<3; i++) get_pointer_current_fp_ampere(lev, i)->FillBoundary(Geom(lev).periodicity());
}

void WarpX::HybridPICSolveE (DtType a_dt_type)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        HybridPICSolveE(lev, a_dt_type);
    }
    FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);
}

void WarpX::HybridPICSolveE (const int lev, DtType a_dt_type)
{
    WARPX_PROFILE("WarpX::HybridPICSolveE()");

    HybridPICSolveE(lev, PatchType::fine, a_dt_type);
    if (lev > 0)
    {
        amrex::Abort(Utils::TextMsg::Err(
        "HybridPICSolveE: Only one level implemented for hybrid-PIC solver."));
    }
}

void WarpX::HybridPICSolveE (const int lev, PatchType patch_type, DtType a_dt_type)
{
    // Solve E field in regular cells
    // The first half step uses t=n quantities, the second half t=n+1/2
    // quantities and the full step uses t=n+1 quantities
    if (a_dt_type == DtType::FirstHalf) {
        m_fdtd_solver_fp[lev]->HybridPICSolveE(
            Efield_fp[lev], current_fp_ampere[lev], current_fp_temp[lev],
            Bfield_fp[lev], rho_fp_temp[lev], electron_pressure_fp[lev],
            m_edge_lengths[lev], lev, m_hybrid_pic_model, a_dt_type
        );
    }
    else if (a_dt_type == DtType::SecondHalf) {
        m_fdtd_solver_fp[lev]->HybridPICSolveE(
            Efield_fp[lev], current_fp_ampere[lev], current_fp[lev],
            Bfield_fp[lev], rho_fp_temp[lev], electron_pressure_fp[lev],
            m_edge_lengths[lev], lev, m_hybrid_pic_model, a_dt_type
        );
    }
    else {
        m_fdtd_solver_fp[lev]->HybridPICSolveE(
            Efield_fp[lev], current_fp_ampere[lev], current_fp_temp[lev],
            Bfield_fp[lev], rho_fp[lev], electron_pressure_fp[lev],
            m_edge_lengths[lev], lev, m_hybrid_pic_model, a_dt_type
        );
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
}

void WarpX::CalculateElectronPressure(DtType a_dt_type)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        CalculateElectronPressure(lev, a_dt_type);
    }
}

void WarpX::CalculateElectronPressure(const int lev, DtType a_dt_type)
{
    // The full step uses rho^{n+1}, otherwise use the old or averaged
    // charge density.
    if (a_dt_type == DtType::Full) {
        m_hybrid_pic_model->FillElectronPressureMF(
            electron_pressure_fp[lev], rho_fp[lev]
        );
    } else {
        m_hybrid_pic_model->FillElectronPressureMF(
            electron_pressure_fp[lev], rho_fp_temp[lev]
        );
    }
    ApplyElectronPressureBoundary(lev, PatchType::fine);
    electron_pressure_fp[lev]->FillBoundary(Geom(lev).periodicity());
}
