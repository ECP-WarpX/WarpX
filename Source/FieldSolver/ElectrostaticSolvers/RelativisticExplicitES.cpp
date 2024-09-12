#include "WarpX.H"

#include "RelativisticExplicitES.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"

using namespace amrex;

void RelativisticExplicitES::InitData () {
    auto & warpx = WarpX::GetInstance();
    bool prepare_field_solve = (WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic);
    // check if any of the particle containers have initialize_self_fields = True
    for (auto const& species : warpx.GetPartContainer()) {
        prepare_field_solve |= species->initialize_self_fields;
    }
    prepare_field_solve |= m_poisson_boundary_handler->m_boundary_potential_specified;

    if (prepare_field_solve) {
        m_poisson_boundary_handler->DefinePhiBCs(warpx.Geom(0));
    }
}

void RelativisticExplicitES::ComputeSpaceChargeField (
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_fp,
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_cp,
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& charge_buf,
    amrex::Vector< std::unique_ptr<amrex::MultiFab> >& phi_fp,
    MultiParticleContainer& mpc,
    MultiFluidContainer* mfl,
    amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Efield_fp,
    amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Bfield_fp
) {
    WARPX_PROFILE("RelativisticExplicitES::ComputeSpaceChargeField");
    amrex::ignore_unused(rho_fp, rho_cp, phi_fp, mfl);

    const bool always_run_solve = (WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic);

    // Loop over the species and add their space-charge contribution to E and B.
    // Note that the fields calculated here does not include the E field
    // due to simulation boundary potentials
    for (auto const& species : mpc) {
        if (always_run_solve || (species->initialize_self_fields)) {
            AddSpaceChargeField(charge_buf, *species, Efield_fp, Bfield_fp);
        }
    }

    // Add the field due to the boundary potentials
    if (always_run_solve || (m_poisson_boundary_handler->m_boundary_potential_specified))
    {
        AddBoundaryField(Efield_fp);
    }
}

void RelativisticExplicitES::AddSpaceChargeField (
    amrex::Vector<std::unique_ptr<amrex::MultiFab> >& charge_buf,
    WarpXParticleContainer& pc,
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield_fp,
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield_fp)
{
    WARPX_PROFILE("RelativisticExplicitES::AddSpaceChargeField");

    if (pc.getCharge() == 0) { return; }

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::n_rz_azimuthal_modes == 1,
                                     "Error: RZ electrostatic only implemented for a single mode");
#endif

    auto & warpx = WarpX::GetInstance();

    // Allocate fields for charge and potential
    Vector<std::unique_ptr<MultiFab> > rho(num_levels);
    Vector<std::unique_ptr<MultiFab> > rho_coarse(num_levels); // Used in order to interpolate between levels
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    // Use number of guard cells used for local deposition of rho
    const amrex::IntVect ng = warpx.get_ng_depos_rho();
    for (int lev = 0; lev < num_levels; lev++) {
        BoxArray nba = warpx.boxArray(lev);
        nba.surroundingNodes();
        rho[lev] = std::make_unique<MultiFab>(nba, warpx.DistributionMap(lev), 1, ng);
        rho[lev]->setVal(0.);
        phi[lev] = std::make_unique<MultiFab>(nba, warpx.DistributionMap(lev), 1, 1);
        phi[lev]->setVal(0.);
        if (lev > 0) {
            // For MR levels: allocated the coarsened version of rho
            BoxArray cba = nba;
            cba.coarsen(warpx.refRatio(lev-1));
            rho_coarse[lev] = std::make_unique<MultiFab>(cba, warpx.DistributionMap(lev), 1, ng);
            rho_coarse[lev]->setVal(0.);
            if (charge_buf[lev]) {
                charge_buf[lev]->setVal(0.);
            }
        }
    }
    // Deposit particle charge density (source of Poisson solver)
    // The options below are identical to those in MultiParticleContainer::DepositCharge
    bool const local = true;
    bool const reset = false;
    bool const apply_boundary_and_scale_volume = true;
    bool const interpolate_across_levels = false;
    if ( !pc.do_not_deposit) {
        pc.DepositCharge(rho, local, reset, apply_boundary_and_scale_volume,
                              interpolate_across_levels);
    }
    warpx.SyncRho(rho, rho_coarse, charge_buf); // Apply filter, perform MPI exchange, interpolate across levels

    // Get the particle beta vector
    bool const local_average = false; // Average across all MPI ranks
    std::array<ParticleReal, 3> beta_pr = pc.meanParticleVelocity(local_average);
    std::array<Real, 3> beta;
    for (int i=0 ; i < static_cast<int>(beta.size()) ; i++) {
        beta[i] = beta_pr[i]/PhysConst::c; // Normalize
    }

    // Compute the potential phi, by solving the Poisson equation
    computePhi( rho, phi, beta, pc.self_fields_required_precision,
                pc.self_fields_absolute_tolerance, pc.self_fields_max_iters,
                pc.self_fields_verbosity );

    // Compute the corresponding electric and magnetic field, from the potential phi
    computeE( Efield_fp, phi, beta );
    computeB( Bfield_fp, phi, beta );

}

void RelativisticExplicitES::AddBoundaryField (amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield_fp)
{
    WARPX_PROFILE("RelativisticExplicitES::AddBoundaryField");

    auto & warpx = WarpX::GetInstance();

    // Allocate fields for charge and potential
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rho(num_levels);
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > phi(num_levels);
    // Use number of guard cells used for local deposition of rho
    const amrex::IntVect ng = warpx.get_ng_depos_rho();
    for (int lev = 0; lev < num_levels; lev++) {
        BoxArray nba = warpx.boxArray(lev);
        nba.surroundingNodes();
        rho[lev] = std::make_unique<amrex::MultiFab>(nba, warpx.DistributionMap(lev), 1, ng);
        rho[lev]->setVal(0.);
        phi[lev] = std::make_unique<amrex::MultiFab>(nba, warpx.DistributionMap(lev), 1, 1);
        phi[lev]->setVal(0.);
    }

    // Set the boundary potentials appropriately
    setPhiBC(phi, warpx.gett_new(0));

    // beta is zero for boundaries
    const std::array<Real, 3> beta = {0._rt};

    // Compute the potential phi, by solving the Poisson equation
    computePhi( rho, phi, beta, self_fields_required_precision,
                self_fields_absolute_tolerance, self_fields_max_iters,
                self_fields_verbosity );

    // Compute the corresponding electric field, from the potential phi.
    computeE( Efield_fp, phi, beta );
}
