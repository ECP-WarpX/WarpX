/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Evolve/WarpXDtType.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"
#include "Fluids/MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <ablastr/fields/MultiFabRegister.H>

using namespace amrex;

void WarpX::HybridPICEvolveFields ()
{
    using ablastr::fields::Direction;
    using ablastr::fields::va2vm;

    WARPX_PROFILE("WarpX::HybridPICEvolveFields()");

    // The below deposition is hard coded for a single level simulation
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        finest_level == 0,
        "Ohm's law E-solve only works with a single level.");

    // The particles have now been pushed to their t_{n+1} positions.
    // Perform charge deposition in component 0 of rho_fp at t_{n+1}.
    mypc->DepositCharge(m_fields.get_mr_levels("rho_fp", finest_level), 0._rt);
    // Perform current deposition at t_{n+1/2}.
    mypc->DepositCurrent(m_fields.get_mr_levels_alldirs("current_fp", finest_level), dt[0], -0.5_rt * dt[0]);

    // Deposit cold-relativistic fluid charge and current
    if (do_fluid_species) {
        int const lev = 0;
        myfl->DepositCharge(m_fields, *m_fields.get("rho_fp", lev), lev);
        myfl->DepositCurrent(m_fields,
            *m_fields.get("current_fp", Direction{0}, lev),
            *m_fields.get("current_fp", Direction{1}, lev),
            *m_fields.get("current_fp", Direction{2}, lev),
            lev);
    }

    // Synchronize J and rho:
    // filter (if used), exchange guard cells, interpolate across MR levels
    // and apply boundary conditions
    SyncCurrentAndRho();

    // SyncCurrent does not include a call to FillBoundary, but it is needed
    // for the hybrid-PIC solver since current values are interpolated to
    // a nodal grid
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            m_fields.get("current_fp", Direction{idim}, lev)->FillBoundary(Geom(lev).periodicity());
        }
    }

    // Get requested number of substeps to use
    const int sub_steps = m_hybrid_pic_model->m_substeps;

    // Get the external current
    m_hybrid_pic_model->GetCurrentExternal(
        m_fields.get_mr_levels_alldirs("edge_lengths", finest_level));

    // Reference hybrid-PIC multifabs
    auto& rho_fp_temp = m_hybrid_pic_model->rho_fp_temp;
    auto& current_fp_temp = m_hybrid_pic_model->current_fp_temp;

    // During the above deposition the charge and current density were updated
    // so that, at this time, we have rho^{n} in rho_fp_temp, rho{n+1} in the
    // 0'th index of `rho_fp`, J_i^{n-1/2} in `current_fp_temp` and J_i^{n+1/2}
    // in `current_fp`.

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
                0.5_rt, *m_fields.get("current_fp", Direction{idim}, lev), 0,
                0, 1, current_fp_temp[lev][idim]->nGrowVect()
            );
        }
    }

    ablastr::fields::MultiLevelVectorField const& Bfield_fp = m_fields.get_mr_levels_alldirs("Bfield_fp", finest_level);

    // Push the B field from t=n to t=n+1/2 using the current and density
    // at t=n, while updating the E field along with B using the electron
    // momentum equation
    for (int sub_step = 0; sub_step < sub_steps; sub_step++)
    {
        m_hybrid_pic_model->BfieldEvolveRK(
            m_fields.get_mr_levels_alldirs("Bfield_fp", finest_level),
            m_fields.get_mr_levels_alldirs("Efield_fp", finest_level),
            va2vm(current_fp_temp), amrex::GetVecOfPtrs(rho_fp_temp),
            m_fields.get_mr_levels_alldirs("edge_lenghts", finest_level),
            0.5_rt/sub_steps*dt[0],
            DtType::FirstHalf, guard_cells.ng_FieldSolver,
            WarpX::sync_nodal_points
        );
    }

    // Average rho^{n} and rho^{n+1} to get rho^{n+1/2} in rho_fp_temp
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Perform a linear combination of values in the 0'th index (1 comp)
        // of rho^{n} and rho^{n+1} (with 0.5 prefactors), writing
        // the result into the 0'th index of `rho_fp_temp[lev]`
        MultiFab::LinComb(
            *rho_fp_temp[lev], 0.5_rt, *rho_fp_temp[lev], 0,
            0.5_rt, *m_fields.get("rho_fp", lev), 0, 0, 1, rho_fp_temp[lev]->nGrowVect()
        );
    }

    // Now push the B field from t=n+1/2 to t=n+1 using the n+1/2 quantities
    for (int sub_step = 0; sub_step < sub_steps; sub_step++)
    {
        m_hybrid_pic_model->BfieldEvolveRK(
            m_fields.get_mr_levels_alldirs("Bfield_fp", finest_level),
            m_fields.get_mr_levels_alldirs("Efield_fp", finest_level),
            m_fields.get_mr_levels_alldirs("current_fp", finest_level),
            amrex::GetVecOfPtrs(rho_fp_temp),
            m_fields.get_mr_levels_alldirs("edge_lenghts", finest_level),
            0.5_rt/sub_steps*dt[0],
            DtType::SecondHalf, guard_cells.ng_FieldSolver,
            WarpX::sync_nodal_points
        );
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
                2._rt, *m_fields.get("current_fp", Direction{idim}, lev), 0,
                0, 1, current_fp_temp[lev][idim]->nGrowVect()
            );
        }
    }

    // Calculate the electron pressure at t=n+1
    m_hybrid_pic_model->CalculateElectronPressure();

    // Update the E field to t=n+1 using the extrapolated J_i^n+1 value
    m_hybrid_pic_model->CalculateCurrentAmpere(
        m_fields.get_mr_levels_alldirs("Bfield_fp", finest_level),
        m_fields.get_mr_levels_alldirs("edge_lengths", finest_level));
    m_hybrid_pic_model->HybridPICSolveE(
        m_fields.get_mr_levels_alldirs("Efield_fp", finest_level),
        va2vm(current_fp_temp),
        m_fields.get_mr_levels_alldirs("Bfield_fp", finest_level),
        m_fields.get_mr_levels("rho_fp", finest_level),
        m_fields.get_mr_levels_alldirs("edge_lengths", finest_level), false
    );
    FillBoundaryE(guard_cells.ng_FieldSolver, WarpX::sync_nodal_points);

    // Copy the rho^{n+1} values to rho_fp_temp and the J_i^{n+1/2} values to
    // current_fp_temp since at the next step those values will be needed as
    // rho^{n} and J_i^{n-1/2}.
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // copy 1 component value starting at index 0 to index 0
        MultiFab::Copy(*rho_fp_temp[lev], *m_fields.get("rho_fp", lev),
                        0, 0, 1, rho_fp_temp[lev]->nGrowVect());
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::Copy(*current_fp_temp[lev][idim], *m_fields.get("current_fp", Direction{idim}, lev),
                           0, 0, 1, current_fp_temp[lev][idim]->nGrowVect());
        }
    }
}

void WarpX::HybridPICDepositInitialRhoAndJ ()
{
    using ablastr::fields::va2vm;

    auto& rho_fp_temp = m_hybrid_pic_model->rho_fp_temp;
    auto& current_fp_temp = m_hybrid_pic_model->current_fp_temp;
    mypc->DepositCharge(amrex::GetVecOfPtrs(rho_fp_temp), 0._rt);
    mypc->DepositCurrent(va2vm(current_fp_temp), dt[0], 0._rt);
    SyncRho(amrex::GetVecOfPtrs(rho_fp_temp), m_fields.get_mr_levels("rho_cp", finest_level), m_fields.get_mr_levels("rho_buf", finest_level));
    SyncCurrent("current_fp");
    for (int lev=0; lev <= finest_level; ++lev) {
        // SyncCurrent does not include a call to FillBoundary, but it is needed
        // for the hybrid-PIC solver since current values are interpolated to
        // a nodal grid
        current_fp_temp[lev][0]->FillBoundary(Geom(lev).periodicity());
        current_fp_temp[lev][1]->FillBoundary(Geom(lev).periodicity());
        current_fp_temp[lev][2]->FillBoundary(Geom(lev).periodicity());

        ApplyRhofieldBoundary(lev, rho_fp_temp[lev].get(), PatchType::fine);
        // Set current density at PEC boundaries, if needed.
        ApplyJfieldBoundary(
            lev, current_fp_temp[lev][0].get(),
            current_fp_temp[lev][1].get(),
            current_fp_temp[lev][2].get(),
            PatchType::fine
        );
    }
}
