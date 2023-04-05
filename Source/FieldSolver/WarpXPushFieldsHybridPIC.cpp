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

    // Perform charge deposition in component 0 of rho_fp
    mypc->DepositCharge(rho_fp, 0._rt);
    // Perform current deposition
    mypc->DepositCurrent(current_fp, dt[0], 0._rt);

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

    // TODO: insert Runge-Kutta integration logic for B update instead
    // of the substep update used here - can test with small timestep using
    // this simpler implementation

    // Note: E^{n} is recalculated with the accurate J_i^{n} since at the end
    // of the last step we had to "guess" J_i^{n}. It also needs to be
    // recalculated to include the resistivity before evolving B.

    // Firstly J_i^{n} is calculated as the average of J_i^{n-1/2} and J_i^{n+1/2}
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*current_fp[lev][0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Extract field data for this grid/tile
            Array4<Real const> const& Jx_adv = current_fp[lev][0]->const_array(mfi);
            Array4<Real const> const& Jy_adv = current_fp[lev][1]->const_array(mfi);
            Array4<Real const> const& Jz_adv = current_fp[lev][2]->const_array(mfi);
            Array4<Real> const& Jx = current_fp_temp[lev][0]->array(mfi);
            Array4<Real> const& Jy = current_fp_temp[lev][1]->array(mfi);
            Array4<Real> const& Jz = current_fp_temp[lev][2]->array(mfi);

            // Extract tileboxes for which to loop
            Box const& tjx  = mfi.tilebox(current_fp_temp[lev][0]->ixType().toIntVect());
            Box const& tjy  = mfi.tilebox(current_fp_temp[lev][1]->ixType().toIntVect());
            Box const& tjz  = mfi.tilebox(current_fp_temp[lev][2]->ixType().toIntVect());

            amrex::ParallelFor(tjx, tjy, tjz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Jx(i, j, k) = 0.5 * (Jx(i, j, k) + Jx_adv(i, j, k));
                },

                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Jy(i, j, k) = 0.5 * (Jy(i, j, k) + Jy_adv(i, j, k));
                },

                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Jz(i, j, k) = 0.5 * (Jz(i, j, k) + Jz_adv(i, j, k));
                }
            );
        }

        // fill ghost cells with appropriate values
        for (int idim = 0; idim < 3; ++idim) {
            current_fp_temp[lev][idim]->FillBoundary(Geom(lev).periodicity());
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

    // Average rho^{n} and rho^{n+1} to get rho^{n+1/2}
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*rho_fp[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            // Extract field data for this grid/tile
            Array4<Real const> const& rho_adv = rho_fp[lev]->const_array(mfi);
            Array4<Real> const& rho = rho_fp_temp[lev]->array(mfi);

            // Extract tilebox for which to loop
            Box const& tb  = mfi.tilebox(rho_fp_temp[lev]->ixType().toIntVect());

            amrex::ParallelFor(tb,
                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    rho(i, j, k) = 0.5 * (rho(i, j, k) + rho_adv(i, j, k));
                }
            );
        }
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

    // Calculate the electron pressure at t=n+1
    CalculateElectronPressure(DtType::Full);

    // Extrapolate the ion current density to t=n+1 to calculate a projected
    // E at t=n+1
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*current_fp[lev][0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Extract field data for this grid/tile
            Array4<Real const> const& Jx_adv = current_fp[lev][0]->const_array(mfi);
            Array4<Real const> const& Jy_adv = current_fp[lev][1]->const_array(mfi);
            Array4<Real const> const& Jz_adv = current_fp[lev][2]->const_array(mfi);
            Array4<Real> const& Jx = current_fp_temp[lev][0]->array(mfi);
            Array4<Real> const& Jy = current_fp_temp[lev][1]->array(mfi);
            Array4<Real> const& Jz = current_fp_temp[lev][2]->array(mfi);

            // Extract tileboxes for which to loop
            Box const& tjx  = mfi.tilebox(current_fp_temp[lev][0]->ixType().toIntVect());
            Box const& tjy  = mfi.tilebox(current_fp_temp[lev][1]->ixType().toIntVect());
            Box const& tjz  = mfi.tilebox(current_fp_temp[lev][2]->ixType().toIntVect());

            amrex::ParallelFor(tjx, tjy, tjz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Jx(i, j, k) = 2.0 * Jx_adv(i, j, k) - Jx(i, j, k);
                },

                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Jy(i, j, k) = 2.0 * Jy_adv(i, j, k) - Jy(i, j, k);
                },

                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Jz(i, j, k) = 2.0 * Jz_adv(i, j, k) - Jz(i, j, k);
                }
            );
        }

        // fill ghost cells with appropriate values
        for (int idim = 0; idim < 3; ++idim) {
            current_fp_temp[lev][idim]->FillBoundary(Geom(lev).periodicity());
        }
    }

    // Update the E field to t=n+1 using the extrapolated J_i^n+1 value
    CalculateCurrentAmpere();
    HybridPICSolveE(DtType::Full);

    // Copy the rho^{n+1} values to rho_fp_temp and the J_i^{n+1/2} values to
    // current_fp_temp since at the next step those values will be needed as
    // rho^{n} and J_i^{n-1/2}.
    for (int lev = 0; lev <= finest_level; ++lev)
    {
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

void WarpX::CalculateCurrentAmpere (int lev)
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

void WarpX::HybridPICSolveE (int lev, DtType a_dt_type)
{
    WARPX_PROFILE("WarpX::HybridPICSolveE()");

    HybridPICSolveE(lev, PatchType::fine, a_dt_type);
    if (lev > 0)
    {
        amrex::Abort(Utils::TextMsg::Err(
        "HybridPICSolveE: Only one level implemented for hybrid-PIC solver."));
        HybridPICSolveE(lev, PatchType::coarse, a_dt_type);
    }
}

void WarpX::HybridPICSolveE (int lev, PatchType patch_type, DtType a_dt_type)
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
