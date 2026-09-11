/* Copyright 2026 Prabhat Kumar
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Fields.H"
#include "ThetaImplicitHybrid.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "EmbeddedBoundary/Enabled.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/ElectronPressureFlux.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/ExternalVectorPotential.H"
#include "Particles/MultiParticleContainer.H"
#include "WarpX.H"
#include <ablastr/utils/Communication.H>
#include <ablastr/coarsen/sample.H>

#include <string>

using warpx::fields::FieldType;
using namespace amrex::literals;

void ThetaImplicitHybrid::Define (WarpX* const a_WarpX, bool /*from_restart*/)
{
    BL_PROFILE("ThetaImplicitHybrid::Define()");

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "ThetaImplicitHybrid object is already defined!");

    m_WarpX = a_WarpX;
    m_num_amr_levels = 1;
    // E0 for the mass-matrix Jacobian is saved inside ComputeRHS (see the comment
    // there); the default SaveE in PreLinearSolve would save the wrong field.
    m_scheme_saves_E0 = true;

    m_hybrid_pic_model = m_WarpX->get_pointer_HybridPICModel();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_hybrid_pic_model != nullptr,
        "ThetaImplicitHybrid solver requires hybrid PIC model to be defined");
    // With the electron energy equation on, the implicit scheme handles the
    // transport, compression and Joule heating through the in-loop pe advance
    // (discretely energy-paired); only the symmetric Q_ei ion-electron exchange
    // is applied once per step (it must kick the ion particles). The
    // include_joule_heating flag is therefore inert implicitly, and the
    // Joule-redirect-to-ions option is not supported.
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !(m_hybrid_pic_model->m_solve_electron_energy_equation &&
          m_hybrid_pic_model->m_joule_redirect_to_ions),
        "The Joule redirect-to-ions option is not supported with the "
        "theta-implicit hybrid solver (Joule heating enters through the "
        "in-loop energy pairing).");

    /// Set flag for external fields from vector potentials
    m_add_external_fields = m_hybrid_pic_model->m_add_external_fields;

    {
        const amrex::ParmParse pp_impl("implicit_evolve");
        pp_impl.query("pe_newton_unknown", m_pe_unknown);
        // enthalpy-flux discretization of the in-loop pe advance: the van Albada MUSCL
        // face flux needs the collocated (nodal J) Cartesian grid and is the default
        // there; central differences elsewhere
#if defined(WARPX_DIM_RZ)
        const bool pe_adv_muscl_ok = false;
#else
        const bool pe_adv_muscl_ok = m_WarpX->m_fields.get(
            FieldType::current_fp, ablastr::fields::Direction{2}, 0)->ixType().nodeCentered();
#endif
        std::string pe_adv_name = pe_adv_muscl_ok ? "vanalbada" : "central";
        pp_impl.query("pe_advection", pe_adv_name);
        if (pe_adv_name == "central") { m_pe_advection = 0; }
        else if (pe_adv_name == "vanalbada") { m_pe_advection = 1; }
        else {
            WARPX_ABORT_WITH_MESSAGE(
                "implicit_evolve.pe_advection = " + pe_adv_name +
                " is not valid; options: central, vanalbada");
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_pe_advection == 0 || pe_adv_muscl_ok,
            "implicit_evolve.pe_advection = vanalbada requires the collocated Cartesian grid");
        pp_impl.query("filter_push_fields", m_filter_push_fields);
        if (m_filter_push_fields) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::use_filter,
                "implicit_evolve.filter_push_fields requires warpx.use_filter = 1");
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_add_external_fields,
                "implicit_evolve.filter_push_fields: external fields are not "
                "supported (the external-field work ledger is unfiltered)");
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_WarpX->Geom(0).isPeriodic(d),
                    "implicit_evolve.filter_push_fields requires a fully "
                    "periodic domain (binomial-filter self-adjointness at "
                    "walls is not handled)");
            }
        }
        pp_impl.query("pe_ue_cap_fac", m_pe_ue_cap_fac);
        if (m_pe_unknown) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_hybrid_pic_model->m_solve_electron_energy_equation &&
                !m_hybrid_pic_model->m_implicit_use_algebraic_closure,
                "implicit_evolve.pe_newton_unknown requires the in-loop electron "
                "energy equation (not the algebraic closure)");
            // pe-row scale: E_pe = -grad(pe)/(e n) ~ pe/(e n0 dx), so 1/(q_e n0_ref dx_min)
            // makes the pe residual commensurate with E in the solver norms
            amrex::Real dx_min = m_WarpX->Geom(0).CellSize(0);
            for (int d = 1; d < AMREX_SPACEDIM; ++d) {
                dx_min = std::min(dx_min, m_WarpX->Geom(0).CellSize(d));
            }
            m_pe_scale = 1.0_rt /
                (PhysConst::q_e * m_hybrid_pic_model->m_n0_ref * dx_min);
        }
        // the preconditioner reads the pe-row scale through the hybrid model
        m_hybrid_pic_model->m_pe_newton_scale = m_pe_unknown ? m_pe_scale : 1.0_rt;
    }
    if (m_pe_unknown) {
        m_E.Define( m_WarpX, "Efield_fp", "hybrid_electron_pressure_fp" );
    } else {
        m_E.Define( m_WarpX, "Efield_fp" );
    }
    m_Eold.Define( m_E );

    // Define B_old MultiFabs
    using ablastr::fields::Direction;
    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        const auto& Bfp_x = m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{0}, lev);
        const auto& dm = Bfp_x->DistributionMap();
        const amrex::IntVect ngb = Bfp_x->nGrowVect();

        for (int dir = 0; dir < 3; ++dir) {
            const auto& ba = m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{dir}, lev)->boxArray();
            m_WarpX->m_fields.alloc_init(FieldType::B_old, Direction{dir}, lev, ba, dm, 1, ngb, 0.0_rt);
        }
    }

    const amrex::ParmParse pp("implicit_evolve");
    pp.query("theta", m_theta);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_theta >= 0.5 && m_theta <= 1.0,
        "theta parameter must be between 0.5 and 1.0");

    parseNonlinearSolverParams( pp );

    if (m_use_mass_matrices) {
        // lagged mass matrices (see DepositMassMatricesThisIter)
        pp.query("mass_matrices_deposit_interval", m_mass_matrices_deposit_interval);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_mass_matrices_deposit_interval >= 0,
            "implicit_evolve.mass_matrices_deposit_interval must be >= 0");
        pp.query("mass_matrices_step_interval", m_mass_matrices_step_interval);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_mass_matrices_step_interval >= 1,
            "implicit_evolve.mass_matrices_step_interval must be >= 1");
    }

    m_nlsolver->Define(m_E, this);

    if (m_pe_unknown) {
        // the pe row must pass through the preconditioner (identity at
        // minimum) or the right-preconditioned operator is singular in it
        const PreconditionerType pc_type = m_nlsolver->GetPreconditionerType();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            pc_type == PreconditionerType::none ||
            pc_type == PreconditionerType::pc_hybrid_pic,
            "implicit_evolve.pe_newton_unknown: only pc_type none or "
            "pc_hybrid_pic pass the pe row through the preconditioner");
    }

    if (m_use_mass_matrices) { InitializeMassMatrices(); }

    m_is_defined = true;
}

void ThetaImplicitHybrid::PrintParameters () const
{
    BL_PROFILE("ThetaImplicitHybrid::PrintParameters()");

    if (!m_WarpX->Verbose()) { return; }
    amrex::Print() << "\n";
    amrex::Print() << "-----------------------------------------------------------\n";
    amrex::Print() << "-------- THETA IMPLICIT HYBRID PIC SOLVER PARAMETERS ------\n";
    amrex::Print() << "-----------------------------------------------------------\n";
    amrex::Print() << "Time-bias parameter theta:           " << m_theta << "\n";
    if (m_use_mass_matrices) {
        amrex::Print() << "mass matrices deposit interval:      "
                       << m_mass_matrices_deposit_interval << "\n";
        amrex::Print() << "mass matrices step interval:         "
                       << m_mass_matrices_step_interval << "\n";
    }
    PrintBaseImplicitSolverParameters();
    m_nlsolver->PrintParams();
    amrex::Print() << "-----------------------------------------------------------\n\n";
}

int ThetaImplicitHybrid::OneStep ( const amrex::Real  start_time,
                                    const amrex::Real  a_dt,
                                    const int          a_step )
{
    BL_PROFILE("ThetaImplicitHybrid::OneStep()");

    m_dt = a_dt;

    // Handle external field splitting: work with internal fields during the solve
    if (m_add_external_fields) {
        m_hybrid_pic_model->m_external_vector_potential->UpdateHybridExternalFields(
            start_time, 0.5_rt * m_dt);
        SubtractExternalEfield();
        SubtractExternalBfield();
    }

    // Save particle state at t^n
    m_WarpX->SaveParticlesAtImplicitStepStart();

    // Save E^n (and pe^n when pe is a Newton unknown)
    if (m_pe_unknown) {
        amrex::MultiFab* pe =
            m_WarpX->m_fields.get(FieldType::hybrid_electron_pressure_fp, 0);
        if (!m_pe_old) {
            // first use: the pressure field holds its (seeded) initialization
            m_pe_old = std::make_unique<amrex::MultiFab>(
                pe->boxArray(), pe->DistributionMap(), pe->nComp(), pe->nGrowVect());
            amrex::MultiFab::Copy(*m_pe_old, *pe, 0, 0, pe->nComp(), pe->nGrowVect());
            m_pe_old->FillBoundary(m_WarpX->Geom(0).periodicity());
        }
        m_Eold.Copy(FieldType::Efield_fp,
                    FieldType::hybrid_electron_pressure_fp);
        m_Eold.getScalarVec()[0]->mult(m_pe_scale, 0, 1);
    } else {
        m_Eold.Copy(FieldType::Efield_fp);
    }

    // Save B^n
    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        const ablastr::fields::VectorField Bfp = m_WarpX->m_fields.get_alldirs(FieldType::Bfield_fp, lev);
        ablastr::fields::VectorField B_old = m_WarpX->m_fields.get_alldirs(FieldType::B_old, lev);
        for (int n = 0; n < 3; ++n) {
            amrex::MultiFab::Copy(*B_old[n], *Bfp[n], 0, 0,
                                  B_old[n]->nComp(), B_old[n]->nGrowVect());
        }
    }

    // Initial guess: E^{n+θ} = E^n
    m_E.Copy(m_Eold);

    // Lagged mass matrices across steps (see DepositMassMatricesThisIter):
    // deposit on the first step of the run/restart and on every k-th step.
    if (m_use_mass_matrices) {
        m_mm_deposit_this_step = (m_mass_matrices_step_interval <= 1) || !m_mm_deposited_once
                                 || (a_step % m_mass_matrices_step_interval == 0);
    }

    // Solve nonlinear system for E^{n+θ} (and eventually Pe^{n+θ})
    m_nlsolver->Solve( m_E, m_Eold, start_time, m_dt, a_step );

    if (m_mm_deposit_this_step) { m_mm_deposited_once = true; }

    const int exit_status = m_nlsolver->GetExitStatus();
    if (exit_status < 0) { return exit_status; }

    // Update WarpX fields to t^{n+θ}
    UpdateWarpXFields( m_E, start_time );
    if (m_pe_unknown) {
        // the accepted Newton iterate is pe^{n+theta}; snapshot it for the t^{n+1}
        // extrapolation in FinishFieldUpdate (the in-loop path snapshots in the advance)
        AMREX_ALWAYS_ASSERT(m_pe_theta);
        const amrex::MultiFab* pe =
            m_WarpX->m_fields.get(FieldType::hybrid_electron_pressure_fp, 0);
        amrex::MultiFab::Copy(*m_pe_theta, *pe, 0, 0,
                              pe->nComp(), pe->nGrowVect());
    }
    m_WarpX->reduced_diags->ComputeDiagsMidStep(a_step);

    // Advance particles from t^{n+1/2} to t^{n+1}
    m_WarpX->FinishImplicitParticleUpdate( start_time + m_dt );

    // Advance fields (including the electron pressure state) from t^{n+θ} to t^{n+1}
    FinishFieldUpdate( start_time + m_dt );

    // Electron energy equation: the transport/compression/Joule part of the
    // update already happened inside the Newton solve (in-loop pe advance,
    // energy-paired). What remains is the symmetric Q_ei ion-electron
    // collisional exchange, which must kick the ion particles and therefore
    // runs once per step here, on the t^{n+1} state.
    if (m_hybrid_pic_model->m_solve_electron_energy_equation &&
        m_hybrid_pic_model->m_include_temperature_relaxation) {
        // Without the Q_ei relaxation there is nothing left to do: transport,
        // compression and Joule heating are already handled by the in-loop pe
        // advance, so the energy equation reduces exactly to the gamma-law
        // closure path (no end-of-step redistribute/deposits needed).
        //
        // FinishImplicitParticleUpdate moved particles to t^{n+1} but did not
        // Redistribute them into their valid cells. QDSMCApplyIonHeating does a
        // per-ion NGP lookup into a zero-guard coefficient MultiFab, so ions left
        // in guard cells would read out of bounds. Redistribute first, matching
        // the explicit path's precondition (it redistributes before the QDSMC step).
        m_WarpX->GetPartContainer().Redistribute();
        m_WarpX->GetPartContainer().DepositCharge(
            m_WarpX->m_fields.get_mr_levels(FieldType::rho_fp, m_num_amr_levels - 1),
            0._rt);
        // Sync rho^{n+1} (DepositCharge does not): fold guard-cell deposits into
        // the valid (incl. periodic) nodes so n_e is unbiased at the boundaries.
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            amrex::MultiFab* rf = m_WarpX->m_fields.get(FieldType::rho_fp, lev);
            ablastr::utils::communication::SumBoundary(
                *rf, 0, rf->nComp(), rf->nGrowVect(), rf->nGrowVect(),
                WarpX::do_single_precision_comms, m_WarpX->Geom(lev).periodicity());
            ablastr::utils::communication::FillBoundary(
                *rf, rf->nGrowVect(), WarpX::do_single_precision_comms,
                m_WarpX->Geom(lev).periodicity(), true);
        }
        // Per-species charge densities rho_fp_<spec> at t^{n+1}: the QDSMC Joule
        // and Q_ei sources read them for the species fractions
        // f_s = rho_s / Sigma_t rho_t. The explicit path deposits them every step
        // in HybridPICDepositRhoAndJ; the implicit path does not call that, so
        // without this deposit they stay frozen at their initialization values
        // (stale f_s, and f_s = 0 in cells the plasma has since moved into).
        // Deposited unscaled in RZ (apply_boundary_and_scale_volume = false) to
        // match the explicit convention -- the 2*pi*r factors cancel in f_s.
        {
            auto & mypc = m_WarpX->GetPartContainer();
            for (auto const & spec : mypc.GetSpeciesNames()) {
                auto & pc = mypc.GetParticleContainerFromName(spec);
                if (pc.getCharge() == 0._prt) { continue; }
                pc.DepositCharge(
                    m_WarpX->m_fields.get_mr_levels("rho_fp_" + spec, m_num_amr_levels - 1),
                    /*local*/false, /*reset*/true,
                    /*apply_boundary_and_scale_volume*/false,
                    /*interpolate_across_levels*/false);
            }
        }
        // Deposit the per-species ion temperature T_<nm> for the Q_ei relaxation.
        // The explicit path fills it in HybridPICDepositRhoAndJ; the implicit path
        // does not call that, so it must deposit here (on the redistributed t^{n+1}
        // particles) or the Q_ei exchange reads a stale T_i.
        m_WarpX->GetPartContainer().DepositTemperatures(m_WarpX->m_fields, 0._rt);
        // In-loop integration: transport, compression and Joule heating were
        // advanced inside the Newton solve (energy-paired); only the symmetric
        // Q_ei ion-electron exchange remains, applied on T_e^{n+1} synced from
        // the in-loop pe^{n+1} and the freshly deposited rho^{n+1}.
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            m_hybrid_pic_model->FillTeFromPe(lev);
            m_hybrid_pic_model->ApplyIonElectronEnergyExchange(lev, m_dt);
            m_hybrid_pic_model->FillPeFromTe(lev);
        }
        // Roll the in-loop pressure state so the next step starts from the
        // relaxed pe^{n+1}.
        if (m_pe_old) {
            amrex::MultiFab* pe =
                m_WarpX->m_fields.get(FieldType::hybrid_electron_pressure_fp, 0);
            amrex::MultiFab::Copy(*m_pe_old, *pe, 0, 0, pe->nComp(), pe->nGrowVect());
            m_pe_old->FillBoundary(m_WarpX->Geom(0).periodicity());
        }
    } else {
        // No Q_ei step: still mirror T_e = P_e/(n_e k_B) from the in-loop
        // pe^{n+1} so the "Te" diagnostic tracks the evolving pressure --
        // it is otherwise only filled at initialization and would dump as a
        // stale uniform value. Diagnostic-only (rho_fp here is the last
        // solver-state deposit, an O(theta dt) old density).
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            m_hybrid_pic_model->FillTeFromPe(lev);
        }
    }

    return exit_status;
}

void ThetaImplicitHybrid::ComputeRHS ( WarpXSolverVec&        a_RHS,
                                       const WarpXSolverVec&  a_E,
                                       amrex::Real            start_time,
                                       int                    a_nl_iter,
                                       bool                   a_from_jacobian )
{
    BL_PROFILE("ThetaImplicitHybrid::ComputeRHS()");

    UpdateWarpXFields( a_E, start_time );

    const amrex::Real theta_time = start_time + m_theta * m_dt;

    ablastr::fields::MultiLevelVectorField Efield_fp =
        m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::Efield_fp, m_num_amr_levels - 1);
    ablastr::fields::MultiLevelVectorField Bfield_fp =
        m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::Bfield_fp, m_num_amr_levels - 1);
    ablastr::fields::MultiLevelVectorField current_fp =
        m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::current_fp, m_num_amr_levels - 1);
    ablastr::fields::MultiLevelScalarField rho_fp =
        m_WarpX->m_fields.get_mr_levels(FieldType::rho_fp, m_num_amr_levels - 1);

    m_hybrid_pic_model->CalculatePlasmaCurrent(Bfield_fp, m_WarpX->GetEBUpdateEFlag());

    if (m_use_mass_matrices_jacobian && a_from_jacobian && m_Ji_save[0]) {
        // Use the ion current frozen at the last nonlinear evaluation, so the push
        // field is a pure function of the Newton variable (see m_Ji_save).
        for (int n = 0; n < 3; ++n) {
            amrex::MultiFab::Copy(*current_fp[0][n], *m_Ji_save[n], 0, 0,
                                  m_Ji_save[n]->nComp(), m_Ji_save[n]->nGrowVect());
        }
    }

    // Particles are pushed with the Newton iterate itself, minus the dissipative
    // part of Ohm's law: E* = a_E - D, D = eta*J_p - eta_h*nabla^2(J_p)
    // (Stanier et al. JCP 2019, Eq. (1); E* = a_E for eta = eta_h = 0). Pushing with the iterate gives the residual a true
    // Jacobian through the particle response -- in particular the electrostatic
    // limit (B = 0) is degenerate with any recomputed push field, which would
    // not depend on the solver variable at all.
    SubtractDissipativeEFromPushField();

    m_WarpX->ApplyFillBoundaryE();

    if (m_add_external_fields) {
        m_hybrid_pic_model->m_external_vector_potential->UpdateHybridExternalFields(
            theta_time, 0.5_rt * m_dt);
        AddExternalBfield();
        AddExternalEfield();
    }

    if (!a_from_jacobian && m_use_mass_matrices_jacobian) {
        // Save the push field E0 for the mass-matrix linear model J = J0 + MM*(E - E0).
        // This must be the same field the linear stage sees at this point of the
        // evaluation (the resistivity-free Ohm's-law E, incl. external fields), not the
        // full Ohm's-law E that Efield_fp holds after ComputeRHS: saving the latter
        // (the default SaveE in PreLinearSolve) puts an O(||R||) offset into MM*(E-E0)
        // and stalls Newton once the fluctuation amplitude grows.
        SaveE();
        if (WarpX::use_filter) {
            // PreRHSOp filters Efield_fp in place before the linear stage contracts the
            // mass matrices with it, so E0 must be the filtered push field as well
            m_WarpX->ApplyFilterMF(
                m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::Efield_fp_save, 0), 0);
        }
        // Also capture the ion current that produced this push field (see m_Ji_save).
        for (int n = 0; n < 3; ++n) {
            const amrex::MultiFab& J = *current_fp[0][n];
            if (!m_Ji_save[n]) {
                m_Ji_save[n] = std::make_unique<amrex::MultiFab>(
                    J.boxArray(), J.DistributionMap(), J.nComp(), J.nGrowVect());
            }
            amrex::MultiFab::Copy(*m_Ji_save[n], J, 0, 0, J.nComp(), J.nGrowVect());
        }
    }

    if (m_filter_push_fields) {
        // conservative smoothing: the particles gather the filtered push field, everything
        // downstream keeps the unfiltered registry (Jacobian probes included)
        FilterPushFieldsSwap(true);
    }
    PreRHSOp( theta_time, a_nl_iter, a_from_jacobian );
    if (m_filter_push_fields) {
        FilterPushFieldsSwap(false);
    }

    {
        // Make the Ohm's-law rho a pure function of the Newton iterate: component 0
        // (deposited at entry positions) carries the previous evaluation's particle
        // state, which pollutes the finite-difference Jacobian matvec at O(hysteresis/eps)
        // in problems where E is deposition-noise dominated (e.g. the electrostatic
        // limit). The post-push component is deposited fresh from x^n each evaluation
        // (Stanier et al. use the half-time moments in Ohm's law for the same reason).
        // At the nonlinear fixed point the two components coincide.
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            amrex::MultiFab* rho = m_WarpX->m_fields.get(FieldType::rho_fp, lev);
            const int nc = rho->nComp()/2;
            const int c_new = rho->nComp() - nc;
            // NOTE: not MultiFab::Copy(*rho, *rho, ...) -- self-copy is UB in AMReX
            for (amrex::MFIter mfi(*rho); mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.growntilebox();
                auto const& a = rho->array(mfi);
                amrex::ParallelFor(bx, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    a(i,j,k,n) = a(i,j,k,c_new + n);
                });
            }
            // The deposit+SumBoundary leaves rho ghosts stale; the electron-pressure
            // advance reads rho through ghosts at box faces, so refresh them here.
            rho->FillBoundary(m_WarpX->Geom(lev).periodicity());
        }
    }

    if (m_use_mass_matrices_jacobian) {
        if (!a_from_jacobian) {
            // Nonlinear evaluation: J and rho are now scaled and synced; capture the
            // linearization base for the rho response (see m_J_base in the header).
            CaptureJRhoBase();
        } else if (m_J_base[0]) {
            // Linear stage: rho = rho_base - (dt/2) div(J - J_base). This replaces both
            // the frozen-rho and retained-rho semantics of the MM linear stage.
            ApplyRhoResponseFromDivJ();
        }
    }

    if (m_add_external_fields) {
        SubtractExternalBfield();
        SubtractExternalEfield();
    }

    // --- Compute full Ohm's law E for Faraday update ---
    // The gamma-law electron pressure is advanced to pe^{n+theta} inside every
    // residual evaluation (in-loop), with its work terms discretely paired to
    // the ion push and Faraday ledgers -- the implicit scheme's
    // implementation of the gamma-law closure with a discrete electron
    // energy ledger. With implicit_use_algebraic_closure on, the pressure
    // is instead re-evaluated algebraically from the iterate's rho (the
    // exact counterpart of the explicit path's default closure): no pe
    // state, no work pairing, and none of the pairing's floored-edge
    // artifacts (Yee Hall-work residual, Cartesian m=4 separatrix layer).
    if (m_hybrid_pic_model->m_implicit_use_algebraic_closure) {
        m_hybrid_pic_model->CalculateElectronPressure();
    } else {
        AdvanceElectronPressure( a_from_jacobian, theta_time );
    }

    m_hybrid_pic_model->HybridPICSolveE(
        Efield_fp, current_fp, Bfield_fp, rho_fp,
        m_WarpX->GetEBUpdateEFlag(),
        true, true   // with resistivity and ∇Pe included (∇Pe is curl-free so doesn't affect Faraday, but needed for self-consistent Newton residual)
    );

    // EB: the masked Ohm solve leaves the Newton iterate at covered locations (a zero,
    // singular residual row). Stamp E = 0 there, the explicit path's convention, so the
    // covered rows become identity rows F = E of the Jacobian.
    if (EB::enabled()) {
        using warpx::fields::FieldType;
        auto const& eb_flags = m_WarpX->GetEBUpdateEFlag()[0];
        for (int dim = 0; dim < 3; ++dim) {
            amrex::MultiFab* Ed = m_WarpX->m_fields.get(
                FieldType::Efield_fp, ablastr::fields::Direction{dim}, 0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(*Ed, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.tilebox();
                auto const& e = Ed->array(mfi);
                auto const& f = eb_flags[dim]->const_array(mfi);
                amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        if (f(i,j,k) == 0) { e(i,j,k) = 0.0_rt; }
                    });
            }
        }
    }

    m_WarpX->ApplyFillBoundaryE();

    // RHS = E_ohm - E_old
    if (m_pe_unknown) {
        a_RHS.Copy(FieldType::Efield_fp, FieldType::hybrid_electron_pressure_fp);
        // pe row: the pressure field holds the iterate (read by Ohm's law above); the
        // residual uses the one-evaluation theta update m_pe_rhs, RHS_pe = kappa*(pe_rhs - pe^n)
        amrex::MultiFab::Copy(*a_RHS.getScalarVec()[0], *m_pe_rhs, 0, 0, 1,
                              amrex::IntVect::TheZeroVector());
        a_RHS.getScalarVec()[0]->mult(m_pe_scale, 0, 1);
    } else {
        a_RHS.Copy(FieldType::Efield_fp);
    }
    a_RHS.linComb(1.0, a_RHS, -1.0, m_Eold);
}

void ThetaImplicitHybrid::CaptureJRhoBase ()
{
    using warpx::fields::FieldType;
    using ablastr::fields::Direction;

    const int lev = 0;
    const ablastr::fields::VectorField J = m_WarpX->m_fields.get_alldirs(FieldType::current_fp, lev);
    const amrex::MultiFab* rho = m_WarpX->m_fields.get(FieldType::rho_fp, lev);

    for (int n = 0; n < 3; ++n) {
        if (!m_J_base[n]) {
            m_J_base[n] = std::make_unique<amrex::MultiFab>(
                J[n]->boxArray(), J[n]->DistributionMap(), J[n]->nComp(), J[n]->nGrowVect());
        }
        amrex::MultiFab::Copy(*m_J_base[n], *J[n], 0, 0, J[n]->nComp(), J[n]->nGrowVect());
    }
    if (!m_rho_base) {
        m_rho_base = std::make_unique<amrex::MultiFab>(
            rho->boxArray(), rho->DistributionMap(), rho->nComp(), rho->nGrowVect());
    }
    amrex::MultiFab::Copy(*m_rho_base, *rho, 0, 0, rho->nComp(), rho->nGrowVect());
}

void ThetaImplicitHybrid::ApplyRhoResponseFromDivJ ()
{
    // rho = rho_base - (dt/2) * div(J - J_base), applied to the component of rho that
    // the Ohm's-law solve reads (component 0). J components are interpolated to nodes
    // (rho is nodal in the hybrid model) and differenced centrally; in RZ the m=0
    // cylindrical divergence is used, with the axis limit (1/r)d(r Jr)/dr -> 2 dJr/dr.
    // Higher azimuthal-mode components are left at their base values.
    using namespace amrex::literals;
    using warpx::fields::FieldType;
    using ablastr::fields::Direction;
    using namespace ablastr::coarsen::sample;

    const int lev = 0;
    const ablastr::fields::VectorField J = m_WarpX->m_fields.get_alldirs(FieldType::current_fp, lev);
    amrex::MultiFab* rho = m_WarpX->m_fields.get(FieldType::rho_fp, lev);

    // Start from the base rho (all components)
    amrex::MultiFab::Copy(*rho, *m_rho_base, 0, 0, rho->nComp(), rho->nGrowVect());

    const amrex::Geometry& geom = m_WarpX->Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    [[maybe_unused]] const amrex::Real rmin = geom.ProbLo(0);
    [[maybe_unused]] const amrex::Real dr = geom.CellSize(0);

    const amrex::GpuArray<int, 3> Jx_stag = m_hybrid_pic_model->Jx_IndexType;
    const amrex::GpuArray<int, 3> Jy_stag = m_hybrid_pic_model->Jy_IndexType;
    const amrex::GpuArray<int, 3> Jz_stag = m_hybrid_pic_model->Jz_IndexType;
    const amrex::GpuArray<int, 3> nodal   = {1, 1, 1};
    const amrex::GpuArray<int, 3> coarsen = {1, 1, 1};
    amrex::ignore_unused(Jy_stag);

    const amrex::Real half_dt = 0.5_rt * m_dt;

    // The response is physical only where plasma exists; in near-floor (vacuum) cells
    // the update would be noise that the 1/n factors of Ohm's law amplify. Restrict the
    // update to cells safely above the density floor and clamp the result at the floor.
    const amrex::Real rho_floor = m_hybrid_pic_model->m_n_floor * PhysConst::q_e;

    // Nodal domain box; the central div stencil (with nodal interpolation of J) reads
    // one cell beyond each node, which exceeds the J guard cells at domain-edge nodes.
    // Interior nodes use central differences; the RZ axis uses a one-sided radial
    // difference; remaining edge nodes keep the base rho (zero response there).
    const amrex::Box dom_nodal = amrex::convert(geom.Domain(), amrex::IntVect::TheNodeVector());
    const amrex::Dim3 dlo = amrex::lbound(dom_nodal);
    const amrex::Dim3 dhi = amrex::ubound(dom_nodal);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*rho, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::Array4<amrex::Real>       const& rho_arr = rho->array(mfi);
        amrex::Array4<amrex::Real const> const& Jx  = J[0]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jz  = J[2]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jx0 = m_J_base[0]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jz0 = m_J_base[2]->const_array(mfi);
#if defined(WARPX_DIM_3D)
        amrex::Array4<amrex::Real const> const& Jy  = J[1]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jy0 = m_J_base[1]->const_array(mfi);
#endif

        const amrex::Box tb = mfi.tilebox(amrex::IntVect::TheNodeVector());

        amrex::ParallelFor(tb, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            // nodal current response at neighboring nodes
            auto dJx_n = [&] (int ii, int jj, int kk) {
                return Interp(Jx, Jx_stag, nodal, coarsen, ii, jj, kk, 0)
                     - Interp(Jx0, Jx_stag, nodal, coarsen, ii, jj, kk, 0);
            };
            auto dJz_n = [&] (int ii, int jj, int kk) {
                return Interp(Jz, Jz_stag, nodal, coarsen, ii, jj, kk, 0)
                     - Interp(Jz0, Jz_stag, nodal, coarsen, ii, jj, kk, 0);
            };

            amrex::Real div = 0._rt;
#if defined(WARPX_DIM_RZ)
            if (j <= dlo.y || j >= dhi.y || i >= dhi.x) { return; }
            if (i > dlo.x) {
                const amrex::Real r = rmin + i*dr;
                div += (dJx_n(i+1,j,k) - dJx_n(i-1,j,k)) * 0.5_rt * dxi[0]
                     + dJx_n(i,j,k)/r;
            } else {
                // axis: (1/r) d(r Jr)/dr -> 2 dJr/dr, one-sided
                div += 2._rt * (dJx_n(i+1,j,k) - dJx_n(i,j,k)) * dxi[0];
            }
            div += (dJz_n(i,j+1,k) - dJz_n(i,j-1,k)) * 0.5_rt * dxi[1];
#elif defined(WARPX_DIM_XZ)
            if (i <= dlo.x || i >= dhi.x || j <= dlo.y || j >= dhi.y) { return; }
            div += (dJx_n(i+1,j,k) - dJx_n(i-1,j,k)) * 0.5_rt * dxi[0];
            div += (dJz_n(i,j+1,k) - dJz_n(i,j-1,k)) * 0.5_rt * dxi[1];
#elif defined(WARPX_DIM_1D_Z)
            if (i <= dlo.x || i >= dhi.x) { return; }
            div += (dJz_n(i+1,j,k) - dJz_n(i-1,j,k)) * 0.5_rt * dxi[0];
#elif defined(WARPX_DIM_3D)
            auto dJy_n = [&] (int ii, int jj, int kk) {
                return Interp(Jy, Jy_stag, nodal, coarsen, ii, jj, kk, 0)
                     - Interp(Jy0, Jy_stag, nodal, coarsen, ii, jj, kk, 0);
            };
            if (i <= dlo.x || i >= dhi.x || j <= dlo.y || j >= dhi.y ||
                k <= dlo.z || k >= dhi.z) { return; }
            div += (dJx_n(i+1,j,k) - dJx_n(i-1,j,k)) * 0.5_rt * dxi[0];
            div += (dJy_n(i,j+1,k) - dJy_n(i,j-1,k)) * 0.5_rt * dxi[1];
            div += (dJz_n(i,j,k+1) - dJz_n(i,j,k-1)) * 0.5_rt * dxi[2];
#else
            amrex::ignore_unused(i, j, k, dJx_n, dJz_n, dxi, half_dt, dlo, dhi);
#endif
            const amrex::Real rho0v = rho_arr(i,j,k,0);
            if (rho0v > 10._rt * rho_floor) {
                rho_arr(i,j,k,0) = amrex::max(rho0v - half_dt * div, rho_floor);
            }
        });
    }

    // refresh rho guards for downstream interpolations
    rho->FillBoundary(geom.periodicity());
}

void ThetaImplicitHybrid::SubtractDissipativeEFromPushField ()
{
    // D = E_Ohm(with dissipation) - E_Ohm(without) evaluated from the same
    // (B^{n+theta}, rho, pe) state: the ideal, Hall and grad-pe parts cancel
    // exactly, leaving eta*J_p - eta_h*nabla^2(J_p) with the identical stencils,
    // interpolations, floors and axis handling as the residual's Ohm solve --
    // by construction, for any resistivity model. The FD-solver-level entry is
    // used so no boundary condition is applied to Efield_fp as a side effect.
    using namespace amrex::literals;
    using warpx::fields::FieldType;

    if (m_hybrid_pic_model->m_eta_expression == "0.0" &&
        !m_hybrid_pic_model->m_include_hyper_resistivity_term) { return; }

    const int lev = 0;
    const ablastr::fields::VectorField E =
        m_WarpX->m_fields.get_alldirs(FieldType::Efield_fp, lev);
    const ablastr::fields::VectorField Ji =
        m_WarpX->m_fields.get_alldirs(FieldType::current_fp, lev);
    const amrex::MultiFab* rho = m_WarpX->m_fields.get(FieldType::rho_fp, lev);
    const amrex::MultiFab* pe =
        m_WarpX->m_fields.get(FieldType::hybrid_electron_pressure_fp, lev);
    const ablastr::fields::VectorField B =
        m_WarpX->m_fields.get_alldirs(FieldType::Bfield_fp, lev);

    for (int n = 0; n < 3; ++n) {
        if (!m_D[n]) {
            m_D[n] = std::make_unique<amrex::MultiFab>(
                E[n]->boxArray(), E[n]->DistributionMap(), E[n]->nComp(), E[n]->nGrowVect());
            m_E_work[n] = std::make_unique<amrex::MultiFab>(
                E[n]->boxArray(), E[n]->DistributionMap(), E[n]->nComp(), E[n]->nGrowVect());
            // The masked Ohm solves below skip EB-covered (flag = 0) locations,
            // which therefore retain their allocation-time content forever:
            // without this initialization, D at covered nodes is uninitialized
            // arena memory that pollutes E* = a_E - D and the Newton residual.
            m_D[n]->setVal(0.0_rt);
            m_E_work[n]->setVal(0.0_rt);
        }
    }
    const ablastr::fields::VectorField D    = {m_D[0].get(), m_D[1].get(), m_D[2].get()};
    const ablastr::fields::VectorField Ework = {m_E_work[0].get(), m_E_work[1].get(), m_E_work[2].get()};

    ablastr::fields::VectorField Jp =
        m_WarpX->m_fields.get_alldirs(FieldType::hybrid_current_fp_plasma, lev);
    auto& eb_update_E = m_WarpX->GetEBUpdateEFlag()[lev];
    auto* fdtd = m_WarpX->get_pointer_fdtd_solver_fp(lev);
    // with dissipation (solve_for_Faraday = true), into D
    fdtd->HybridPICSolveE(D, Jp, Ji, B, *rho, *pe, eb_update_E, lev,
                          m_hybrid_pic_model, true, true);
    // without dissipation, into Ework
    fdtd->HybridPICSolveE(Ework, Jp, Ji, B, *rho, *pe, eb_update_E, lev,
                          m_hybrid_pic_model, false, true);

    for (int n = 0; n < 3; ++n) {
        amrex::MultiFab::Subtract(*m_D[n], *m_E_work[n], 0, 0, m_D[n]->nComp(), 0);
        m_D[n]->FillBoundary(m_WarpX->Geom(lev).periodicity());
        // E* = a_E - D (all components; the push and the pairing read this)
        amrex::MultiFab::Subtract(*E[n], *m_D[n], 0, 0, E[n]->nComp(), 0);
    }
}

void ThetaImplicitHybrid::FilterPushFieldsSwap (const bool a_apply)
{
    // Conservative smoothing: PreRHSOp binomial-filters Efield_fp in place for the gather;
    // save E* before and restore it after so Ohm's law and the pe work pairing keep the
    // unfiltered field (B needs no treatment, v x B does no work).
    using warpx::fields::FieldType;
    const int lev = 0;
    const ablastr::fields::VectorField E =
        m_WarpX->m_fields.get_alldirs(FieldType::Efield_fp, lev);
    for (int n = 0; n < 3; ++n) {
        amrex::MultiFab& Emf = *E[n];
        if (a_apply) {
            if (!m_E_unfiltered[n]) {
                m_E_unfiltered[n] = std::make_unique<amrex::MultiFab>(
                    Emf.boxArray(), Emf.DistributionMap(),
                    Emf.nComp(), Emf.nGrowVect());
            }
            amrex::MultiFab::Copy(*m_E_unfiltered[n], Emf, 0, 0,
                                  Emf.nComp(), Emf.nGrowVect());
        } else {
            amrex::MultiFab::Copy(Emf, *m_E_unfiltered[n], 0, 0,
                                  Emf.nComp(), Emf.nGrowVect());
        }
    }
}

void ThetaImplicitHybrid::AdvanceElectronPressure ( const bool a_from_jacobian,
                                                    const amrex::Real a_theta_time )
{
    // pe^{n+theta} = pe^n - theta*dt * [ div(ue pe^n) + (gamma-1) pe^n div(ue) ],
    // with ue = (J_i - J_net)/(e n) evaluated at nodes from the freshly deposited ion
    // current and the plasma (net) current. Forward evaluation in pe is O(dt^2) for
    // the half step and keeps the update an explicit pure function of the iterate.
    using namespace amrex::literals;
    using warpx::fields::FieldType;
    using ablastr::fields::Direction;
    using namespace ablastr::coarsen::sample;

    const int lev = 0;
    amrex::MultiFab* pe = m_WarpX->m_fields.get(FieldType::hybrid_electron_pressure_fp, lev);
    const amrex::MultiFab* rho = m_WarpX->m_fields.get(FieldType::rho_fp, lev);
    const ablastr::fields::VectorField J =
        m_WarpX->m_fields.get_alldirs(FieldType::current_fp, lev);
    const ablastr::fields::VectorField Jp =
        m_WarpX->m_fields.get_alldirs(FieldType::hybrid_current_fp_plasma, lev);
    // Efield_fp holds the push field E* = a_E - D (external-field contributions
    // have already been subtracted again): the exact field the ions were pushed
    // with. theta-Faraday consumes a_E = E* + D, whose dissipative part D enters
    // the pairing as the -D.J_p heating term below.
    const ablastr::fields::VectorField Efld =
        m_WarpX->m_fields.get_alldirs(FieldType::Efield_fp, lev);
    for (int n = 0; n < 3; ++n) {
        if (!m_D[n]) {  // no dissipation configured: pair against D = 0
            m_D[n] = std::make_unique<amrex::MultiFab>(
                Efld[n]->boxArray(), Efld[n]->DistributionMap(),
                Efld[n]->nComp(), Efld[n]->nGrowVect());
            m_D[n]->setVal(0.0_rt);
        }
    }
    // Honor include_joule_heating (default off), matching the explicit QDSMC
    // path: with it off, the dissipative work D.J_p (Ohmic + hyper-resistive)
    // is NOT deposited into pe -- it leaves the system as the vacuum/hyper
    // dissipation it is (Faraday still drains it from W_B, so the B-field
    // damping physics is unchanged). Without this, eta_H's grid-scale
    // dissipation at the sharp FRC edge current sheets heats the few
    // electrons there to keV within ~100 steps. The reversible E*.J_e
    // transport/compression pairing is unaffected. Pair against D = 0 by
    // pointing the kernel at zeroed fields.
    const bool jheat = m_hybrid_pic_model->m_include_joule_heating;
    if (!jheat && !m_D_zero[0]) {
        for (int n = 0; n < 3; ++n) {
            m_D_zero[n] = std::make_unique<amrex::MultiFab>(
                Efld[n]->boxArray(), Efld[n]->DistributionMap(),
                Efld[n]->nComp(), Efld[n]->nGrowVect());
            m_D_zero[n]->setVal(0.0_rt);
        }
    }

    const amrex::Geometry& geom = m_WarpX->Geom(lev);

    // The nodal update below interpolates the cell-centered J and J_plasma through
    // guard cells (a node on a box face needs the J value owned by the neighbor box),
    // so their ghosts must be current before the kernel runs. Without this, the two
    // boxes sharing a nodal point compute different pe there and the pressure field
    // becomes multivalued at box seams (breaking the discrete energy pairing).
    for (int d = 0; d < 3; ++d) {
        J[d]->FillBoundary(geom.periodicity());
        Jp[d]->FillBoundary(geom.periodicity());
    }

    if (!m_pe_old) {
#if defined(WARPX_DIM_RZ)
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::n_rz_azimuthal_modes == 1,
            "implicit hybrid in-loop pe advance: the RZ work pairing is implemented for m = 0 only");
#endif
        // first use: the pressure field holds its (seeded) initialization
        m_pe_old = std::make_unique<amrex::MultiFab>(
            pe->boxArray(), pe->DistributionMap(), pe->nComp(), pe->nGrowVect());
        amrex::MultiFab::Copy(*m_pe_old, *pe, 0, 0, pe->nComp(), pe->nGrowVect());
        m_pe_old->FillBoundary(geom.periodicity());
    }
    if (!m_pe_theta) {
        m_pe_theta = std::make_unique<amrex::MultiFab>(
            pe->boxArray(), pe->DistributionMap(), pe->nComp(), pe->nGrowVect());
        m_pe_scratch = std::make_unique<amrex::MultiFab>(
            pe->boxArray(), pe->DistributionMap(), pe->nComp(), pe->nGrowVect());
    }
    if (m_pe_unknown && !m_pe_rhs) {
        m_pe_rhs = std::make_unique<amrex::MultiFab>(
            pe->boxArray(), pe->DistributionMap(), pe->nComp(), pe->nGrowVect());
    }

    const auto dxi = geom.InvCellSizeArray();
    [[maybe_unused]] const amrex::Real rmin = geom.ProbLo(0);
    [[maybe_unused]] const amrex::Real dr = geom.CellSize(0);
    const amrex::Real theta_dt = m_theta * m_dt;
    const amrex::Real gamma = m_hybrid_pic_model->m_gamma;
    const amrex::Real rho_floor = m_hybrid_pic_model->m_n_floor * PhysConst::q_e;
    const amrex::Real q_e = PhysConst::q_e;
    // width of the C1 positivity floor (HybridPeFloor) as a fraction of the floored-adiabat
    // pe: a hard max(pe, 0) in the residual is non-differentiable where floored cells ride
    // it and stalls Newton; the hard clamp remains in the post-step halo pin
    constexpr amrex::Real pe_eps_fac = 0.01_rt;
    const amrex::Real pe_eps = pe_eps_fac * m_hybrid_pic_model->m_n_floor
        * m_hybrid_pic_model->m_elec_temp
        * std::pow(m_hybrid_pic_model->m_n_floor / m_hybrid_pic_model->m_n0_ref,
                   gamma - 1._rt);
    // marker-CFL cap on u_e: the theta-centered fixed point contracts only where
    // |u_e| k_max theta dt < 1, which the fictitious floored-edge u_e = J/(e n_floor)
    // violates (the physical interior u_e sits well below the cap)
    amrex::Real dx_min = geom.CellSize(0);
    for (int d = 1; d < AMREX_SPACEDIM; ++d) {
        dx_min = std::min(dx_min, geom.CellSize(d));
    }
    const amrex::Real ue_cap = m_pe_ue_cap_fac * dx_min / theta_dt;

    const amrex::GpuArray<int, 3> Jx_stag = m_hybrid_pic_model->Jx_IndexType;
    const amrex::GpuArray<int, 3> Jy_stag = m_hybrid_pic_model->Jy_IndexType;
    const amrex::GpuArray<int, 3> Jz_stag = m_hybrid_pic_model->Jz_IndexType;
    const amrex::GpuArray<int, 3> nodal   = {1, 1, 1};
    const amrex::GpuArray<int, 3> coarsen = {1, 1, 1};
    amrex::ignore_unused(Jy_stag);

    const amrex::Box dom_nodal = amrex::convert(geom.Domain(), amrex::IntVect::TheNodeVector());
    const amrex::Dim3 dlo = amrex::lbound(dom_nodal);
    const amrex::Dim3 dhi = amrex::ubound(dom_nodal);
    amrex::GpuArray<bool, 3> is_per = {true, true, true};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { is_per[d] = geom.isPeriodic(d); }

    // Collocated grid: J, E and rho are all nodal and HybridPICSolveE builds the
    // pressure field with CartesianNodalAlgorithm (centered difference at the node,
    // nodal rho) instead of the Yee edge construction. The work pairing below must
    // use the identical stencil or the electron side subtracts a different discrete
    // work than the ions receive (J-correlated leak).
    const bool J_nodal = J[2]->ixType().nodeCentered();

    // Joule heating on the collocated grid: deposit the positive-definite Q (HybridPeJouleQ,
    // see m_Q_diss) instead of -D.Jp, whose hyper-resistive part is sign-indefinite pointwise
    const bool q_posdef = jheat && J_nodal;
    if (q_posdef) {
        if (!m_Q_diss) {
            m_Q_diss = std::make_unique<amrex::MultiFab>(
                pe->boxArray(), pe->DistributionMap(), 1,
                amrex::IntVect::TheZeroVector());
            m_Q_diss->setVal(0.0_rt);
        }
        const ablastr::fields::VectorField Bf =
            m_WarpX->m_fields.get_alldirs(FieldType::Bfield_fp, lev);
        const auto eta_ex  = m_hybrid_pic_model->m_eta;
        const auto etah_ex = m_hybrid_pic_model->m_eta_h;
        const bool inc_hyp = m_hybrid_pic_model->m_include_hyper_resistivity_term;
        const amrex::Real t_now = a_theta_time;
        const amrex::Dim3 qlo = dlo, qhi = dhi;
        const amrex::GpuArray<bool, 3> qper = is_per;
        const auto dxiq = geom.InvCellSizeArray();
        for (amrex::MFIter mfi(*m_Q_diss, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const amrex::Box tb = mfi.tilebox();
            auto const& q   = m_Q_diss->array(mfi);
            auto const& jpx = Jp[0]->const_array(mfi);
            auto const& jpy = Jp[1]->const_array(mfi);
            auto const& jpz = Jp[2]->const_array(mfi);
            auto const& bxa = Bf[0]->const_array(mfi);
            auto const& bya = Bf[1]->const_array(mfi);
            auto const& bza = Bf[2]->const_array(mfi);
            auto const& rr  = rho->const_array(mfi);
            // shared with the explicit fluid pe solver (ElectronPressureFlux.H)
            amrex::ParallelFor(tb,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    q(i,j,k) = HybridPeJouleQ(i, j, k, jpx, jpy, jpz, bxa, bya, bza, rr,
                                              eta_ex, etah_ex, inc_hyp, t_now, dxiq,
                                              qlo, qhi, qper);
                });
        }
    }

    // electron conduction (numeric kappa_e or a kappa_e(rho,Te) expression); pe and rho are
    // nodal on every grid, RZ m = 0 uses the area-weighted radial divergence below
    const bool has_kappa = m_hybrid_pic_model->m_has_kappa_e;
    const bool has_kexpr = m_hybrid_pic_model->m_has_kappa_e_expression;
    const amrex::Real kappa_e = m_hybrid_pic_model->m_kappa_e;
    const auto kappa_ex = m_hybrid_pic_model->m_kappa;
    if (has_kappa) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            pe->nGrowVect().min() >= 1 && rho->nGrowVect().min() >= 1,
            "hybrid_pic_model.kappa_e (electron conduction) needs >= 1 "
            "ghost cell on pe and rho");
    }

    const int pe_adv = m_pe_advection;
    if (pe_adv != 0) {
        // the MUSCL face reconstruction reads the second node beyond a periodic box face
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(pe->nGrowVect().min() >= 2,
            "implicit_evolve.pe_advection = vanalbada needs >= 2 pe ghost cells");
    }
    // face reconstruction parameters, shared with the explicit fluid pe solver
    HybridPeFluxParams pfp;
    pfp.adv = pe_adv;

    // The update must be theta-centered in pe as well: an explicit (forward) pe in the
    // RHS integrates the pressure side of the ion-acoustic oscillation with forward
    // Euler and is numerically unstable (growth ~ exp(omega^2 dt t / 2)). The update is
    // linear in pe, so a short fixed-point iteration (contraction ~ theta*omega*dt)
    // converges the theta-centered value: pe_rhs = (1-theta)*pe^n + theta*pe_iter.
    // Contraction per cycle is ~theta*k_max*Cs*dt; 4 cycles verified sufficient
    // (8 cycles bit-reproduces the energy history on the cold-beam FGI test).
    const int n_pe_iters = m_pe_unknown ? 1 : 4;
    amrex::MultiFab* pe_out = pe;
    if (m_pe_unknown) {
        // pe is a Newton unknown: evaluate the theta update once at the iterate (held by the
        // pressure field for Ohm's law) into m_pe_rhs for the residual row. Non-periodic
        // domain-face nodes keep pe^n, so their rows are identities kappa*(pe - pe^n).
        amrex::MultiFab::Copy(*m_pe_scratch, *pe, 0, 0, pe->nComp(), pe->nGrowVect());
        amrex::MultiFab::Copy(*m_pe_rhs, *m_pe_old, 0, 0, pe->nComp(), pe->nGrowVect());
        // freeze the flux stencil's wall values at pe^n as the in-loop path does: a wall
        // value slaved to the interior iterate feeds the near-wall flux back onto itself
        // (unstable wall layer); the wall rows themselves stay identities
        for (amrex::MFIter mfi(*m_pe_scratch); mfi.isValid(); ++mfi) {
            const amrex::Box vbx = mfi.validbox();
            auto const& s  = m_pe_scratch->array(mfi);
            auto const& p0 = m_pe_old->const_array(mfi);
            const amrex::Dim3 dl = dlo, dh = dhi;
            const amrex::GpuArray<bool, 3> per = is_per;
            amrex::ParallelFor(vbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    amrex::ignore_unused(j, k);
                    const bool face =
                        ((i == dl.x || i == dh.x) && !per[0])
#if (AMREX_SPACEDIM >= 2)
                        || ((j == dl.y || j == dh.y) && !per[1])
#endif
#if (AMREX_SPACEDIM == 3)
                        || ((k == dl.z || k == dh.z) && !per[2])
#endif
                        ;
                    if (face) { s(i,j,k) = p0(i,j,k); }
                });
        }
        pe_out = m_pe_rhs.get();
    } else {
        amrex::MultiFab::Copy(*m_pe_scratch, *m_pe_old, 0, 0, pe->nComp(), pe->nGrowVect());
    }
    for (int pe_it = 0; pe_it < n_pe_iters; ++pe_it) {

    // EB: freeze covered nodes at pe^n (identity rows with pe_unknown, no-op otherwise) so
    // the wall-adjacent Ohm rows do not read an advected floor-density state through
    // grad(pe); the nodal flag is the component-0 E flag of the masked Ohm solves
    const amrex::iMultiFab* eb_pe_flag = EB::enabled()
        ? m_WarpX->GetEBUpdateEFlag()[0][0].get() : nullptr;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*pe, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        amrex::Array4<amrex::Real>       const& pe_arr  = pe_out->array(mfi);
        amrex::Array4<amrex::Real const> const& pe_it_arr = m_pe_scratch->const_array(mfi);
        amrex::Array4<amrex::Real const> const& pe0     = m_pe_old->const_array(mfi);
        amrex::Array4<int const> ebp;
        if (eb_pe_flag) { ebp = eb_pe_flag->const_array(mfi); }
        amrex::Array4<amrex::Real const> const& rho_arr = rho->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jx  = J[0]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jz  = J[2]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jpx = Jp[0]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jpz = Jp[2]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jy  = J[1]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Jpy = Jp[1]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Ex_arr = Efld[0]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Ey_arr = Efld[1]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Ez_arr = Efld[2]->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Dx =
            (jheat ? m_D[0] : m_D_zero[0])->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Dy =
            (jheat ? m_D[1] : m_D_zero[1])->const_array(mfi);
        amrex::Array4<amrex::Real const> const& Dz =
            (jheat ? m_D[2] : m_D_zero[2])->const_array(mfi);
        amrex::Array4<amrex::Real const> Qd;
        if (q_posdef) { Qd = m_Q_diss->const_array(mfi); }

        const amrex::Box tb = mfi.tilebox(amrex::IntVect::TheNodeVector());

        amrex::ParallelFor(tb, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            // EB-covered node: freeze at pe^n (see eb_pe_flag above)
            if (ebp && ebp(i,j,k) == 0) {
                pe_arr(i,j,k) = pe0(i,j,k);
                return;
            }

            // smooth (tanh) marker-CFL cap on u_e, see ue_cap above: a hard min/max makes
            // the residual non-differentiable where floored cells ride the cap and Newton
            // grinds there; tanh keeps the bound and the small-u identity response
            auto uclamp = [&] (amrex::Real u) {
                return ue_cap * std::tanh(u / ue_cap);
            };
            // electron velocity at a node: ue = (J_ion - J_net)/(e n)
            auto ue_x = [&] (int ii, int jj, int kk) {
                const amrex::Real n_ = amrex::max(rho_arr(ii,jj,kk,0), rho_floor);
                return uclamp((Interp(Jx, Jx_stag, nodal, coarsen, ii, jj, kk, 0)
                      - Interp(Jpx, Jx_stag, nodal, coarsen, ii, jj, kk, 0)) / n_);
            };
            auto ue_z = [&] (int ii, int jj, int kk) {
                const amrex::Real n_ = amrex::max(rho_arr(ii,jj,kk,0), rho_floor);
                return uclamp((Interp(Jz, Jz_stag, nodal, coarsen, ii, jj, kk, 0)
                      - Interp(Jpz, Jz_stag, nodal, coarsen, ii, jj, kk, 0)) / n_);
            };
            // The fixed point iterates pe^{n+theta} directly:
            //   pe^{n+theta} = pe^n - theta*dt*RHS(pe^{n+theta}),
            // so the RHS reads the current iterate itself.
            auto pec = [&] (int ii, int jj, int kk) { return pe_it_arr(ii,jj,kk); };
            // flux F = ue * pe_centered and velocity divergence, central differences
            auto Fx = [&] (int ii, int jj, int kk) { return ue_x(ii,jj,kk) * pec(ii,jj,kk); };
            auto Fz = [&] (int ii, int jj, int kk) { return ue_z(ii,jj,kk) * pec(ii,jj,kk); };
            // pe at the face between the p0 and p1 nodes (pm1/p2 = next nodes outward),
            // reconstructed from the upwind side (ElectronPressureFlux.H, shared with the
            // explicit fluid pe solver)
            auto pe_face = [&] (amrex::Real uf, amrex::Real pm1, amrex::Real p0,
                                amrex::Real p1, amrex::Real p2) {
                return HybridPeFace(pfp, uf, pm1, p0, p1, p2);
            };
            // electron heat conduction +div(kappa_e grad(pe/n)) in conservative face-flux
            // form, face conductivity = average of the adjacent nodal values; RZ m = 0 uses
            // the area-weighted radial divergence (axis: r F -> 0)
            auto Tn = [&] (int ii, int jj, int kk) {
                return pec(ii,jj,kk) * q_e
                    / amrex::max(rho_arr(ii,jj,kk,0), rho_floor);
            };
            auto kapn = [&] (int ii, int jj, int kk) {
                if (!has_kexpr) { return kappa_e; }
                return kappa_ex(amrex::max(rho_arr(ii,jj,kk,0), rho_floor),
                                Tn(ii,jj,kk) / q_e);   // Te in eV
            };
            // conduction face flux between two adjacent nodes
            auto Fk = [&] (amrex::Real T0, amrex::Real T1,
                           amrex::Real k0, amrex::Real k1,
                           amrex::Real dxi_f) {
                return 0.5_rt*(k0 + k1) * (T1 - T0) * dxi_f;
            };
            amrex::Real cond = 0.0_rt;
            if (has_kappa) {
#if defined(WARPX_DIM_1D_Z)
                cond = (Fk(Tn(i,j,k),   Tn(i+1,j,k),
                           kapn(i,j,k), kapn(i+1,j,k), dxi[0])
                      - Fk(Tn(i-1,j,k), Tn(i,j,k),
                           kapn(i-1,j,k), kapn(i,j,k), dxi[0])) * dxi[0];
#elif defined(WARPX_DIM_XZ)
                cond = (Fk(Tn(i,j,k),   Tn(i+1,j,k),
                           kapn(i,j,k), kapn(i+1,j,k), dxi[0])
                      - Fk(Tn(i-1,j,k), Tn(i,j,k),
                           kapn(i-1,j,k), kapn(i,j,k), dxi[0])) * dxi[0]
                     + (Fk(Tn(i,j,k),   Tn(i,j+1,k),
                           kapn(i,j,k), kapn(i,j+1,k), dxi[1])
                      - Fk(Tn(i,j-1,k), Tn(i,j,k),
                           kapn(i,j-1,k), kapn(i,j,k), dxi[1])) * dxi[1];
#elif defined(WARPX_DIM_3D)
                cond = (Fk(Tn(i,j,k),   Tn(i+1,j,k),
                           kapn(i,j,k), kapn(i+1,j,k), dxi[0])
                      - Fk(Tn(i-1,j,k), Tn(i,j,k),
                           kapn(i-1,j,k), kapn(i,j,k), dxi[0])) * dxi[0]
                     + (Fk(Tn(i,j,k),   Tn(i,j+1,k),
                           kapn(i,j,k), kapn(i,j+1,k), dxi[1])
                      - Fk(Tn(i,j-1,k), Tn(i,j,k),
                           kapn(i,j-1,k), kapn(i,j,k), dxi[1])) * dxi[1]
                     + (Fk(Tn(i,j,k),   Tn(i,j,k+1),
                           kapn(i,j,k), kapn(i,j,k+1), dxi[2])
                      - Fk(Tn(i,j,k-1), Tn(i,j,k),
                           kapn(i,j,k-1), kapn(i,j,k), dxi[2])) * dxi[2];
#elif defined(WARPX_DIM_RZ)
                const amrex::Real rC = rmin + i*dr;
                const bool onax = (rC < 0.5_rt*dr);
                const amrex::Real Ak = onax ? 0.125_rt*dr*dr : rC*dr;
                const amrex::Real Fkrp = Fk(Tn(i,j,k), Tn(i+1,j,k),
                                            kapn(i,j,k), kapn(i+1,j,k),
                                            dxi[0]);
                const amrex::Real rFkrm = onax ? 0.0_rt :
                    (rC - 0.5_rt*dr) * Fk(Tn(i-1,j,k), Tn(i,j,k),
                                          kapn(i-1,j,k), kapn(i,j,k),
                                          dxi[0]);
                cond = ((rC + 0.5_rt*dr)*Fkrp - rFkrm) / Ak
                     + (Fk(Tn(i,j,k),   Tn(i,j+1,k),
                           kapn(i,j,k), kapn(i,j+1,k), dxi[1])
                      - Fk(Tn(i,j-1,k), Tn(i,j,k),
                           kapn(i,j-1,k), kapn(i,j,k), dxi[1])) * dxi[1];
#endif
            }

            // Cartesian: enthalpy-flux transport plus the work pairing term,
            //   d_t pe = -gamma div(ue pe) + (gamma-1) ue.grad(pe) + (gamma-1) Q.
            //
            // Collocated grid (J_nodal): full-Ohm pairing. The electrons absorb the
            // exact complement of the discrete ion work (+sum E.J_i) and the
            // magnetic-energy change (-sum E.J_amp for theta = 1/2 Faraday with
            // mutually adjoint curls), so W carries -E.J_e = E.(J_i - J_amp) with E
            // the solver iterate that pushed the particles and drives Faraday. This
            // is the continuum (gamma-1)[ue.grad(pe) + Q_Joule] with the conservative
            // Joule form Q = eta J_tot.J_e, and it holds for ANY Ohm's-law contents:
            // the Hall term does no work pointwise at a node, and total energy
            // K_i + U_e + W_B is conserved identically, independent of what E is.
            //
            // Yee grid: legacy pressure-channel-only pairing. E_pe,edge =
            // -UpwardD(pe)/max(rho_edge, floor) is the SAME discrete field
            // HybridPICSolveE builds; each edge's work is split half to each
            // adjacent node so the periodic sum cancels the ion pressure-channel
            // work exactly. The remaining J_net.grad(pe)/(en) term uses centered
            // differences (B-channel pairing not exact on the staggered mesh).
            amrex::Real divF = 0._rt, W = 0._rt;
#if defined(WARPX_DIM_1D_Z)
            if ((i <= dlo.x || i >= dhi.x) && !is_per[0]) { return; }
            if (J_nodal && pe_adv != 0) {
                // limited upwind face flux (conservative: telescopes exactly)
                auto Fz_face = [&] (int m) {
                    const amrex::Real uf = 0.5_rt*(ue_z(m,j,k) + ue_z(m+1,j,k));
                    return uf * pe_face(uf, pec(m-1,j,k), pec(m,j,k),
                                        pec(m+1,j,k), pec(m+2,j,k));
                };
                divF += (Fz_face(i) - Fz_face(i-1)) * dxi[0];
            } else if (J_nodal) {
                divF += (Fz(i+1,j,k) - Fz(i-1,j,k)) * 0.5_rt * dxi[0];
            } else {
                // Yee: conservative edge flux (matches the RZ construction). The
                // central node form would interpolate J through the SECOND ghost
                // ring at box-face nodes, which the deposited current does not
                // have -- the reads land out of bounds and the seam nodes get
                // garbage flux (observed as a box-seam energy leak).
                auto ue_zedge = [&] (int ie) {
                    const amrex::Real rho_e = amrex::max(
                        0.5_rt*(rho_arr(ie,j,k,0) + rho_arr(ie+1,j,k,0)), rho_floor);
                    return uclamp((Jz(ie,j,k) - Jpz(ie,j,k)) / rho_e);
                };
                divF += (ue_zedge(i)   * 0.5_rt*(pe_it_arr(i,j,k)   + pe_it_arr(i+1,j,k))
                       - ue_zedge(i-1) * 0.5_rt*(pe_it_arr(i-1,j,k) + pe_it_arr(i,j,k)))
                      * dxi[0];
            }
            if (J_nodal) {
                // -E*.J_e over all three components (perp components do work at
                // B != 0), minus the dissipative work D.J_p that Faraday drains
                // from W_B (Joule + hyper-resistive heating -> electrons)
                W += Ex_arr(i,j,k) * (Jx(i,j,k) - Jpx(i,j,k))
                   + Ey_arr(i,j,k) * (Jy(i,j,k) - Jpy(i,j,k))
                   + Ez_arr(i,j,k) * (Jz(i,j,k) - Jpz(i,j,k))
                   - (q_posdef ? Qd(i,j,k)
                      : (Dx(i,j,k) * Jpx(i,j,k)
                       + Dy(i,j,k) * Jpy(i,j,k)
                       + Dz(i,j,k) * Jpz(i,j,k)));
            } else {
                // Yee: per-component full-Ohm pairing q_d = E*_d (J_d - Jp_d) - D_d Jp_d
                // at each component's own location; the z-staggered component is
                // half-split to the node (Ex/Ey are nodal in 1D)
                auto qz = [&] (int ie) {
                    return Ez_arr(ie,j,k) * (Jz(ie,j,k) - Jpz(ie,j,k))
                         - Dz(ie,j,k) * Jpz(ie,j,k);
                };
                W += Ex_arr(i,j,k) * (Jx(i,j,k) - Jpx(i,j,k)) - Dx(i,j,k) * Jpx(i,j,k)
                   + Ey_arr(i,j,k) * (Jy(i,j,k) - Jpy(i,j,k)) - Dy(i,j,k) * Jpy(i,j,k)
                   + 0.5_rt * (qz(i-1) + qz(i));
            }
#elif defined(WARPX_DIM_XZ)
            if (((i <= dlo.x || i >= dhi.x) && !is_per[0]) ||
                ((j <= dlo.y || j >= dhi.y) && !is_per[1])) { return; }
            if (J_nodal && pe_adv != 0) {
                // limited upwind face fluxes (conservative: telescope exactly)
                auto Fx_face = [&] (int m) {
                    const amrex::Real uf = 0.5_rt*(ue_x(m,j,k) + ue_x(m+1,j,k));
                    return uf * pe_face(uf, pec(m-1,j,k), pec(m,j,k),
                                        pec(m+1,j,k), pec(m+2,j,k));
                };
                auto Fz_face = [&] (int m) {
                    const amrex::Real uf = 0.5_rt*(ue_z(i,m,k) + ue_z(i,m+1,k));
                    return uf * pe_face(uf, pec(i,m-1,k), pec(i,m,k),
                                        pec(i,m+1,k), pec(i,m+2,k));
                };
                divF += (Fx_face(i) - Fx_face(i-1)) * dxi[0]
                      + (Fz_face(j) - Fz_face(j-1)) * dxi[1];
            } else if (J_nodal) {
                divF += (Fx(i+1,j,k) - Fx(i-1,j,k)) * 0.5_rt * dxi[0]
                      + (Fz(i,j+1,k) - Fz(i,j-1,k)) * 0.5_rt * dxi[1];
            } else {
                // Yee: conservative edge flux; see the 1D branch for why the
                // central node form cannot be used with the deposited current.
                auto ue_xedge = [&] (int ie) {
                    const amrex::Real rho_e = amrex::max(
                        0.5_rt*(rho_arr(ie,j,k,0) + rho_arr(ie+1,j,k,0)), rho_floor);
                    return uclamp((Jx(ie,j,k) - Jpx(ie,j,k)) / rho_e);
                };
                auto ue_zedge = [&] (int je) {
                    const amrex::Real rho_e = amrex::max(
                        0.5_rt*(rho_arr(i,je,k,0) + rho_arr(i,je+1,k,0)), rho_floor);
                    return uclamp((Jz(i,je,k) - Jpz(i,je,k)) / rho_e);
                };
                divF += (ue_xedge(i)   * 0.5_rt*(pe_it_arr(i,j,k)   + pe_it_arr(i+1,j,k))
                       - ue_xedge(i-1) * 0.5_rt*(pe_it_arr(i-1,j,k) + pe_it_arr(i,j,k)))
                      * dxi[0]
                      + (ue_zedge(j)   * 0.5_rt*(pe_it_arr(i,j,k)   + pe_it_arr(i,j+1,k))
                       - ue_zedge(j-1) * 0.5_rt*(pe_it_arr(i,j-1,k) + pe_it_arr(i,j,k)))
                      * dxi[1];
            }
            if (J_nodal) {
                // -E*.J_e over all three components (perp components do work at
                // B != 0), minus the dissipative work D.J_p that Faraday drains
                // from W_B (Joule + hyper-resistive heating -> electrons)
                W += Ex_arr(i,j,k) * (Jx(i,j,k) - Jpx(i,j,k))
                   + Ey_arr(i,j,k) * (Jy(i,j,k) - Jpy(i,j,k))
                   + Ez_arr(i,j,k) * (Jz(i,j,k) - Jpz(i,j,k))
                   - (q_posdef ? Qd(i,j,k)
                      : (Dx(i,j,k) * Jpx(i,j,k)
                       + Dy(i,j,k) * Jpy(i,j,k)
                       + Dz(i,j,k) * Jpz(i,j,k)));
            } else {
                // Yee: per-component full-Ohm pairing, staggered components
                // half-split to the node (Ey is nodal in XZ)
                auto qx = [&] (int ie) {
                    return Ex_arr(ie,j,k) * (Jx(ie,j,k) - Jpx(ie,j,k))
                         - Dx(ie,j,k) * Jpx(ie,j,k);
                };
                auto qz = [&] (int je) {
                    return Ez_arr(i,je,k) * (Jz(i,je,k) - Jpz(i,je,k))
                         - Dz(i,je,k) * Jpz(i,je,k);
                };
                W += 0.5_rt * (qx(i-1) + qx(i))
                   + Ey_arr(i,j,k) * (Jy(i,j,k) - Jpy(i,j,k)) - Dy(i,j,k) * Jpy(i,j,k)
                   + 0.5_rt * (qz(j-1) + qz(j));
            }
#elif defined(WARPX_DIM_RZ)
            // m = 0, Yee-staggered RZ (asserted at first use). The enthalpy flux
            // uses the conservative flux form, so the volume-weighted grid sum
            // telescopes exactly (including through the axis, where [r F] -> 0).
            // W carries the full-Ohm pairing -E*.J_e + eta*Jp^2 evaluated per
            // staggered component at its own location (matching the push-field
            // subtraction), volume-split from the edges to the nodes so that
            // sum_nodes(A W) = sum_edges(A q). The Hall channel and the RZ
            // Faraday ledger pair at interpolation order (not machine-exact,
            // unlike the collocated Cartesian grid).
            amrex::ignore_unused(Fx, Fz, ue_x, ue_z, pe_face);
            if ((j <= dlo.y || j >= dhi.y) && !is_per[1]) { return; }
            if (i >= dhi.x) { return; }                        // outer r boundary
            const amrex::Real r = rmin + i*dr;
            const bool on_axis = (r < 0.5_rt*dr);
            if (i <= dlo.x && !on_axis) { return; }            // annular inner boundary
            const amrex::Real A_i  = on_axis ? 0.125_rt*dr*dr : r*dr;
            const amrex::Real A_ip = (r + 0.5_rt*dr)*dr;       // r-edge at i+1/2
            const amrex::Real A_im = (r - 0.5_rt*dr)*dr;       // r-edge at i-1/2

            // enthalpy flux, conservative form
            auto ue_redge = [&] (int ie) {
                const amrex::Real rho_e = amrex::max(
                    0.5_rt*(rho_arr(ie,j,k,0) + rho_arr(ie+1,j,k,0)), rho_floor);
                return uclamp((Jx(ie,j,k) - Jpx(ie,j,k)) / rho_e);
            };
            auto ue_zedge = [&] (int je) {
                const amrex::Real rho_e = amrex::max(
                    0.5_rt*(rho_arr(i,je,k,0) + rho_arr(i,je+1,k,0)), rho_floor);
                return uclamp((Jz(i,je,k) - Jpz(i,je,k)) / rho_e);
            };
            {
                const amrex::Real Frp =
                    ue_redge(i) * 0.5_rt*(pe_it_arr(i,j,k) + pe_it_arr(i+1,j,k));
                const amrex::Real rFrm = on_axis ? 0.0_rt :
                    (r - 0.5_rt*dr) *
                    ue_redge(i-1) * 0.5_rt*(pe_it_arr(i-1,j,k) + pe_it_arr(i,j,k));
                divF += ((r + 0.5_rt*dr)*Frp - rFrm) / A_i;
                divF += (ue_zedge(j) * 0.5_rt*(pe_it_arr(i,j,k) + pe_it_arr(i,j+1,k))
                       - ue_zedge(j-1) * 0.5_rt*(pe_it_arr(i,j-1,k) + pe_it_arr(i,j,k)))
                      * dxi[1];
            }

            // full-Ohm work pairing per staggered component
            auto q_redge = [&] (int ie) {
                return Ex_arr(ie,j,k) * (Jx(ie,j,k) - Jpx(ie,j,k))
                     - Dx(ie,j,k) * Jpx(ie,j,k);
            };
            auto q_zedge = [&] (int je) {
                return Ez_arr(i,je,k) * (Jz(i,je,k) - Jpz(i,je,k))
                     - Dz(i,je,k) * Jpz(i,je,k);
            };
            {
                const amrex::Real q_theta =
                    Ey_arr(i,j,k) * (Jy(i,j,k) - Jpy(i,j,k)) - Dy(i,j,k) * Jpy(i,j,k);
                const amrex::Real w_r =
                    (0.5_rt * q_redge(i) * A_ip
                     + (on_axis ? 0.0_rt : 0.5_rt * q_redge(i-1) * A_im)) / A_i;
                W += w_r
                   + 0.5_rt * (q_zedge(j-1) + q_zedge(j))
                   + q_theta;
            }
#elif defined(WARPX_DIM_3D)
            auto ue_y = [&] (int ii, int jj, int kk) {
                const amrex::Real n_ = amrex::max(rho_arr(ii,jj,kk,0), rho_floor);
                return uclamp((Interp(Jy, Jy_stag, nodal, coarsen, ii, jj, kk, 0)
                      - Interp(Jpy, Jy_stag, nodal, coarsen, ii, jj, kk, 0)) / n_);
            };
            auto Fy = [&] (int ii, int jj, int kk) { return ue_y(ii,jj,kk) * pec(ii,jj,kk); };
            if (((i <= dlo.x || i >= dhi.x) && !is_per[0]) ||
                ((j <= dlo.y || j >= dhi.y) && !is_per[1]) ||
                ((k <= dlo.z || k >= dhi.z) && !is_per[2])) { return; }
            if (J_nodal && pe_adv != 0) {
                // limited upwind face fluxes (conservative: telescope exactly)
                auto Fx_face = [&] (int m) {
                    const amrex::Real uf = 0.5_rt*(ue_x(m,j,k) + ue_x(m+1,j,k));
                    return uf * pe_face(uf, pec(m-1,j,k), pec(m,j,k),
                                        pec(m+1,j,k), pec(m+2,j,k));
                };
                auto Fy_face = [&] (int m) {
                    const amrex::Real uf = 0.5_rt*(ue_y(i,m,k) + ue_y(i,m+1,k));
                    return uf * pe_face(uf, pec(i,m-1,k), pec(i,m,k),
                                        pec(i,m+1,k), pec(i,m+2,k));
                };
                auto Fz_face = [&] (int m) {
                    const amrex::Real uf = 0.5_rt*(ue_z(i,j,m) + ue_z(i,j,m+1));
                    return uf * pe_face(uf, pec(i,j,m-1), pec(i,j,m),
                                        pec(i,j,m+1), pec(i,j,m+2));
                };
                divF += (Fx_face(i) - Fx_face(i-1)) * dxi[0]
                      + (Fy_face(j) - Fy_face(j-1)) * dxi[1]
                      + (Fz_face(k) - Fz_face(k-1)) * dxi[2];
            } else if (J_nodal) {
                divF += (Fx(i+1,j,k) - Fx(i-1,j,k)) * 0.5_rt * dxi[0]
                      + (Fy(i,j+1,k) - Fy(i,j-1,k)) * 0.5_rt * dxi[1]
                      + (Fz(i,j,k+1) - Fz(i,j,k-1)) * 0.5_rt * dxi[2];
            } else {
                // Yee: conservative edge flux; see the 1D branch for why the
                // central node form cannot be used with the deposited current.
                auto ue_xedge = [&] (int ie) {
                    const amrex::Real rho_e = amrex::max(
                        0.5_rt*(rho_arr(ie,j,k,0) + rho_arr(ie+1,j,k,0)), rho_floor);
                    return uclamp((Jx(ie,j,k) - Jpx(ie,j,k)) / rho_e);
                };
                auto ue_yedge = [&] (int je) {
                    const amrex::Real rho_e = amrex::max(
                        0.5_rt*(rho_arr(i,je,k,0) + rho_arr(i,je+1,k,0)), rho_floor);
                    return uclamp((Jy(i,je,k) - Jpy(i,je,k)) / rho_e);
                };
                auto ue_zedge = [&] (int ke) {
                    const amrex::Real rho_e = amrex::max(
                        0.5_rt*(rho_arr(i,j,ke,0) + rho_arr(i,j,ke+1,0)), rho_floor);
                    return uclamp((Jz(i,j,ke) - Jpz(i,j,ke)) / rho_e);
                };
                divF += (ue_xedge(i)   * 0.5_rt*(pe_it_arr(i,j,k)   + pe_it_arr(i+1,j,k))
                       - ue_xedge(i-1) * 0.5_rt*(pe_it_arr(i-1,j,k) + pe_it_arr(i,j,k)))
                      * dxi[0]
                      + (ue_yedge(j)   * 0.5_rt*(pe_it_arr(i,j,k)   + pe_it_arr(i,j+1,k))
                       - ue_yedge(j-1) * 0.5_rt*(pe_it_arr(i,j-1,k) + pe_it_arr(i,j,k)))
                      * dxi[1]
                      + (ue_zedge(k)   * 0.5_rt*(pe_it_arr(i,j,k)   + pe_it_arr(i,j,k+1))
                       - ue_zedge(k-1) * 0.5_rt*(pe_it_arr(i,j,k-1) + pe_it_arr(i,j,k)))
                      * dxi[2];
            }
            if (J_nodal) {
                // -E*.J_e over all three components (perp components do work at
                // B != 0), minus the dissipative work D.J_p that Faraday drains
                // from W_B (Joule + hyper-resistive heating -> electrons)
                W += Ex_arr(i,j,k) * (Jx(i,j,k) - Jpx(i,j,k))
                   + Ey_arr(i,j,k) * (Jy(i,j,k) - Jpy(i,j,k))
                   + Ez_arr(i,j,k) * (Jz(i,j,k) - Jpz(i,j,k))
                   - (q_posdef ? Qd(i,j,k)
                      : (Dx(i,j,k) * Jpx(i,j,k)
                       + Dy(i,j,k) * Jpy(i,j,k)
                       + Dz(i,j,k) * Jpz(i,j,k)));
            } else {
                // Yee: per-component full-Ohm pairing, each component half-split
                // to the node along its own staggered direction
                auto qx = [&] (int ie) {
                    return Ex_arr(ie,j,k) * (Jx(ie,j,k) - Jpx(ie,j,k))
                         - Dx(ie,j,k) * Jpx(ie,j,k);
                };
                auto qy = [&] (int je) {
                    return Ey_arr(i,je,k) * (Jy(i,je,k) - Jpy(i,je,k))
                         - Dy(i,je,k) * Jpy(i,je,k);
                };
                auto qz = [&] (int ke) {
                    return Ez_arr(i,j,ke) * (Jz(i,j,ke) - Jpz(i,j,ke))
                         - Dz(i,j,ke) * Jpz(i,j,ke);
                };
                W += 0.5_rt * (qx(i-1) + qx(i))
                   + 0.5_rt * (qy(j-1) + qy(j))
                   + 0.5_rt * (qz(k-1) + qz(k));
            }
#else
            amrex::ignore_unused(Fx, Fz, divF, dlo, dhi, is_per, dxi,
                                 pe_face, Tn, kapn, Fk);
#endif
            const amrex::Real pe_new = pe0(i,j,k)
                - theta_dt * (gamma * divF + (gamma - 1._rt) * (W - cond));
            // C1 positivity floor with compact support (exact identity for pe >= 2*eps,
            // hence energy-neutral there), shared with the explicit fluid solver
            pe_arr(i,j,k) = HybridPeFloor(pe_new, pe_eps);
        });
    }

    // Duplicated nodal points on box faces are written redundantly by each box;
    // force them single-valued before they are consumed (Ohm's law, next iterate).
    pe_out->OverrideSync(geom.periodicity());
    pe_out->FillBoundary(geom.periodicity());
    if (pe_it < n_pe_iters - 1) {
        amrex::MultiFab::Copy(*m_pe_scratch, *pe_out, 0, 0, pe->nComp(), pe->nGrowVect());
    }
    } // pe fixed-point iterations

#if defined(WARPX_DIM_RZ)
    {
        // Re-phase the axis-row pressure to the local density through a
        // radially averaged entropy s = pe/rho^gamma (rings 0..2): the axis
        // row's J_e is deposition-noise dominated and the advected pe there
        // dephases from n, closing an anti-restoring feedback loop. Slaving
        // the axis-row pe to n (QDSMC-style rebuild) preserves the entropy
        // evolution while restoring the pressure-density phase lock.
        const amrex::Real gam = m_hybrid_pic_model->m_gamma;
        const amrex::Real rfloor = m_hybrid_pic_model->m_n_floor * PhysConst::q_e;
        const amrex::MultiFab* rho_mf = m_WarpX->m_fields.get(FieldType::rho_fp, lev);
        const amrex::Box dom_n =
            amrex::convert(geom.Domain(), amrex::IntVect::TheNodeVector());
        const int iax = dom_n.smallEnd(0);
        const bool has_axis = (std::abs(geom.ProbLo(0)) < 0.5_rt*geom.CellSize(0));
        if (has_axis) {
            for (amrex::MFIter mfi(*pe_out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const amrex::Box tb = mfi.tilebox(amrex::IntVect::TheNodeVector());
                if (tb.smallEnd(0) > iax) { continue; }
                auto const& pe_a = pe_out->array(mfi);
                auto const& rh   = rho_mf->const_array(mfi);
                const amrex::Box tb0(amrex::IntVect(iax, tb.smallEnd(1)),
                                     amrex::IntVect(iax, tb.bigEnd(1)),
                                     tb.ixType());
                amrex::ParallelFor(tb0, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    amrex::Real sbar = 0._rt;
                    for (int ii = 0; ii <= 2; ++ii) {
                        const amrex::Real nr =
                            amrex::max(rh(i+ii,j,k,0), rfloor);
                        sbar += pe_a(i+ii,j,k) / std::pow(nr, gam);
                    }
                    sbar *= (1._rt/3._rt);
                    const amrex::Real n0 = amrex::max(rh(i,j,k,0), rfloor);
                    pe_a(i,j,k) = sbar * std::pow(n0, gam);
                });
            }
            pe_out->OverrideSync(geom.periodicity());
            pe_out->FillBoundary(geom.periodicity());
        }
    }
#endif

    if (!a_from_jacobian && !m_pe_unknown) {
        amrex::MultiFab::Copy(*m_pe_theta, *pe, 0, 0, pe->nComp(), pe->nGrowVect());
    }
}

void ThetaImplicitHybrid::UpdateWarpXFields ( const WarpXSolverVec&  a_E,
                                                amrex::Real start_time )
{
    BL_PROFILE("ThetaImplicitHybrid::UpdateWarpXFields()");

    const amrex::Real theta_time = start_time + m_theta * m_dt;

    // Set E^{n+θ} in WarpX
    m_WarpX->SetElectricFieldAndApplyBCs( a_E, theta_time );

    // pe as Newton unknown: write the iterate's pe row into the pressure field and fill the
    // domain ghosts (PEC sides included) for Ohm's law. The wall nodes are not rewritten:
    // they are unknowns pinned to pe^n by identity rows (static walls, as in the in-loop path).
    if (m_pe_unknown) {
        amrex::MultiFab* pe =
            m_WarpX->m_fields.get(FieldType::hybrid_electron_pressure_fp, 0);
        amrex::MultiFab::Copy(*pe, *a_E.getScalarVec()[0], 0, 0, 1,
                              amrex::IntVect::TheZeroVector());
        pe->mult(1.0_rt / m_pe_scale, 0, 1);
        pe->FillBoundary(m_WarpX->Geom(0).periodicity());
        m_WarpX->ApplyElectronPressureBoundary(0, PatchType::fine,
                                               /*rewrite_pec_nodes=*/false);
    }

    // Compute B^{n+θ} = B^n - θ·dt·curl(E^{n+θ}) via Faraday's law
    ablastr::fields::MultiLevelVectorField const& B_old =
        m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::B_old, m_num_amr_levels - 1);
    m_WarpX->UpdateMagneticFieldAndApplyBCs( B_old, m_theta * m_dt, start_time );
}


amrex::Array<const amrex::MultiFab*, 3>
ThetaImplicitHybrid::GetBfieldThetaForPC ( const int lev ) const
{
    // During the nonlinear solve, UpdateWarpXFields (called from every
    // residual evaluation) leaves the Bfield_fp registry holding the TOTAL
    // theta-midpoint field B^{n+theta} of the current iterate. Valid only
    // after the first residual evaluation of the current Newton iterate;
    // before that (and between steps) the registry holds the end-of-step
    // totals B^{n+1} (= B^n at the next entry).
    using ablastr::fields::Direction;
    amrex::Array<const amrex::MultiFab*, 3> B = {
        m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{0}, lev),
        m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{1}, lev),
        m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{2}, lev) };
    if (!m_add_external_fields) { return B; }
    // assemble the TOTAL field (see the header doc): Bfield_fp is internal
    // between residual evaluations on this branch
    for (int d = 0; d < 3; ++d) {
        const amrex::MultiFab* Bext = m_WarpX->m_fields.get(
            FieldType::hybrid_B_fp_external, Direction{d}, lev);
        auto& tot = m_B_tot_pc[d];
        if (!tot || tot->boxArray() != B[d]->boxArray()
            || tot->DistributionMap() != B[d]->DistributionMap()) {
            tot = std::make_unique<amrex::MultiFab>(
                B[d]->boxArray(), B[d]->DistributionMap(),
                B[d]->nComp(), B[d]->nGrowVect());
        }
        amrex::MultiFab::Copy(*tot, *B[d], 0, 0,
                              B[d]->nComp(), B[d]->nGrowVect());
        amrex::MultiFab::Add(*tot, *Bext, 0, 0,
                             Bext->nComp(), tot->nGrowVect());
        B[d] = tot.get();
    }
    return B;
}

const amrex::MultiFab*
ThetaImplicitHybrid::GetRhoMidForPC ( const int lev ) const
{
    // The rho_fp registry carries two time slots of WarpX::ncomps
    // components each; consumers read component nComp()/2 (the
    // midpoint-position deposit rho^{n+1/2} of the current iterate,
    // written by every residual evaluation). Valid only after the first
    // residual evaluation of the current Newton iterate.
    return m_WarpX->m_fields.get(FieldType::rho_fp, lev);
}

amrex::Array<const amrex::MultiFab*, 3>
ThetaImplicitHybrid::GetIonCurrentForPC ( const int lev ) const
{
    // The Ohm solve consumes current_fp as the ion (particle) current;
    // it is deposited each residual evaluation and frozen during Jacobian
    // probes, so between preconditioner updates and the GMRES solve it
    // holds exactly the frozen drift-leg coefficient (J - J_i) x delta_B.
    using ablastr::fields::Direction;
    return { m_WarpX->m_fields.get(FieldType::current_fp, Direction{0}, lev),
             m_WarpX->m_fields.get(FieldType::current_fp, Direction{1}, lev),
             m_WarpX->m_fields.get(FieldType::current_fp, Direction{2}, lev) };
}

void ThetaImplicitHybrid::FinishFieldUpdate( amrex::Real end_time )
{
    BL_PROFILE("ThetaImplicitHybrid::FinishFieldUpdate()");

    // B^{n+1}
    ablastr::fields::MultiLevelVectorField const& B_old =
        m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::B_old, 0);
    m_WarpX->FinishMagneticFieldAndApplyBCs( B_old, m_theta, end_time );

    // Add external fields to get total fields at t^{n+1}
    if (m_add_external_fields) {
        m_hybrid_pic_model->m_external_vector_potential->UpdateHybridExternalFields(
            end_time, 0.5_rt * m_dt);
        AddExternalBfield();
        AddExternalEfield();
    }
    // pe^{n+1} = (pe^{n+theta} - (1-theta) pe^n) / theta, then roll the state
    if (m_pe_theta) {
    using namespace amrex::literals;
    using warpx::fields::FieldType;

    const int lev = 0;
    amrex::MultiFab* pe = m_WarpX->m_fields.get(FieldType::hybrid_electron_pressure_fp, lev);
    amrex::MultiFab::LinComb(*pe, 1._rt/m_theta, *m_pe_theta, 0,
                             -(1._rt - m_theta)/m_theta, *m_pe_old, 0,
                             0, pe->nComp(), pe->nGrowVect());
    // Guard against negative overshoot from the extrapolation, and pin
    // below-floor cells to the floored-adiabat constant (the algebraic
    // closure at n_floor -- continuous with the plasma edge and the value
    // the explicit QDSMC's insulating halo holds). Without the pin, the
    // energy-paired work term stores the vacuum-resistivity dissipation
    // eta*J^2 in the halo pe step after step, and T_e = pe/(n_floor k_B)
    // grows without bound. Applied POST-STEP only: pinning inside the
    // Newton residual makes the residual discontinuous in the cells whose
    // re-deposited rho straddles the floor, and Newton stalls.
    {
        const amrex::Real gam = m_hybrid_pic_model->m_gamma;
        const amrex::Real rho_floor = m_hybrid_pic_model->m_n_floor * PhysConst::q_e;
        const amrex::Real pe_vac = m_hybrid_pic_model->m_n_floor
            * m_hybrid_pic_model->m_elec_temp
            * std::pow(m_hybrid_pic_model->m_n_floor / m_hybrid_pic_model->m_n0_ref,
                       gam - 1._rt);
        const amrex::MultiFab* rho = m_WarpX->m_fields.get(FieldType::rho_fp, lev);
        for (amrex::MFIter mfi(*pe); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox();
            auto const& a = pe->array(mfi);
            auto const& r = rho->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                a(i,j,k) = (r(i,j,k,0) <= rho_floor)
                    ? pe_vac : amrex::max(a(i,j,k), 0._rt);
            });
        }
    }
    if (m_pe_unknown) {
        // static walls (see UpdateWarpXFields) held pe^n through the solve; refresh only
        // the domain ghosts for the end-of-step Ohm's-law E below
        pe->FillBoundary(m_WarpX->Geom(lev).periodicity());
        m_WarpX->ApplyElectronPressureBoundary(lev, PatchType::fine,
                                               /*rewrite_pec_nodes=*/false);
    }
    amrex::MultiFab::Copy(*m_pe_old, *pe, 0, 0, pe->nComp(), pe->nGrowVect());
    }

    // E^{n+1}: E is algebraic in the generalized Ohm's law, so evaluate it
    // at the delivered end-of-step state -- total B^{n+1} (externals
    // included above), pe^{n+1}, and the same ion-deposit family the theta
    // stage used. Per-level calls: the multi-level HybridPICSolveE wrapper
    // fires the afterEpush python callback and this must not add a firing.
    {
        using warpx::fields::FieldType;
        if (m_hybrid_pic_model->m_implicit_use_algebraic_closure) {
            m_hybrid_pic_model->CalculateElectronPressure();
        }
        ablastr::fields::MultiLevelVectorField E_fp =
            m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::Efield_fp, m_num_amr_levels - 1);
        ablastr::fields::MultiLevelVectorField J_fp =
            m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::current_fp, m_num_amr_levels - 1);
        ablastr::fields::MultiLevelVectorField B_fp =
            m_WarpX->m_fields.get_mr_levels_alldirs(FieldType::Bfield_fp, m_num_amr_levels - 1);
        ablastr::fields::MultiLevelScalarField r_fp =
            m_WarpX->m_fields.get_mr_levels(FieldType::rho_fp, m_num_amr_levels - 1);
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            m_hybrid_pic_model->HybridPICSolveE(
                E_fp[lev], J_fp[lev], B_fp[lev], *r_fp[lev],
                m_WarpX->GetEBUpdateEFlag()[lev], lev,
                true  /* solve_for_Faraday: include resistivity, eta_H */,
                true  /* keep_grad_pe */);
        }
        // apply the standard E boundary conditions and keep the solver
        // vector consistent with the delivered field
        if (m_pe_unknown) {
            m_E.Copy(FieldType::Efield_fp, FieldType::hybrid_electron_pressure_fp);
            m_E.getScalarVec()[0]->mult(m_pe_scale, 0, 1);
        } else {
            m_E.Copy(FieldType::Efield_fp);
        }
        m_WarpX->SetElectricFieldAndApplyBCs( m_E, end_time );
    }
}

void ThetaImplicitHybrid::AddExternalBfield ()
{
    using ablastr::fields::Direction;

    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            amrex::MultiFab::Add(
                *m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{idim}, lev),
                *m_WarpX->m_fields.get(FieldType::hybrid_B_fp_external, Direction{idim}, lev),
                0, 0, 1,
                m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{idim}, lev)->nGrowVect());
        }
    }
}

void ThetaImplicitHybrid::SubtractExternalBfield ()
{
    using ablastr::fields::Direction;

    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            amrex::MultiFab::Subtract(
                *m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{idim}, lev),
                *m_WarpX->m_fields.get(FieldType::hybrid_B_fp_external, Direction{idim}, lev),
                0, 0, 1,
                m_WarpX->m_fields.get(FieldType::Bfield_fp, Direction{idim}, lev)->nGrowVect());
        }
    }
}

void ThetaImplicitHybrid::AddExternalEfield ()
{
    using ablastr::fields::Direction;

    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            amrex::MultiFab::Add(
                *m_WarpX->m_fields.get(FieldType::Efield_fp, Direction{idim}, lev),
                *m_WarpX->m_fields.get(FieldType::hybrid_E_fp_external, Direction{idim}, lev),
                0, 0, 1,
                m_WarpX->m_fields.get(FieldType::Efield_fp, Direction{idim}, lev)->nGrowVect());
        }
    }
}

void ThetaImplicitHybrid::SubtractExternalEfield ()
{
    using ablastr::fields::Direction;

    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            amrex::MultiFab::Subtract(
                *m_WarpX->m_fields.get(FieldType::Efield_fp, Direction{idim}, lev),
                *m_WarpX->m_fields.get(FieldType::hybrid_E_fp_external, Direction{idim}, lev),
                0, 0, 1,
                m_WarpX->m_fields.get(FieldType::Efield_fp, Direction{idim}, lev)->nGrowVect());
        }
    }
}
