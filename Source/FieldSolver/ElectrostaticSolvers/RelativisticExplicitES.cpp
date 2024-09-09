#include "WarpX.H"

#include "RelativisticExplicitES.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"

void
RelativisticExplicitES::ComputeSpaceChargeField (amrex::Vector<std::unique_ptr<amrex::MultiFab> > charge_buf,
                                                 WarpXParticleContainer& pc, 
                                                 amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield, 
                                                 amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield)
{
    WARPX_PROFILE("WarpX::ComputeSpaceChargeField");

    // Loop over the species and add their space-charge contribution to E and B.
    // Note that the fields calculated here does not include the E field
    // due to simulation boundary potentials
    for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
        WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
        if (species.initialize_self_fields ||
            (electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic)) {
            AddSpaceChargeField(charge_buf, species, Efield, Bfield);
        }
    }

    // Add the field due to the boundary potentials
    if (m_boundary_potential_specified ||
            (electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic)){
        AddBoundaryField(Efield,Bfield);
    }
}

void
RelativisticExplicitES::AddSpaceChargeField (amrex::Vector<std::unique_ptr<amrex::MultiFab> > charge_buf,
                                             WarpXParticleContainer& pc, 
                                             amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield, 
                                             amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield)
{
    WARPX_PROFILE("WarpX::AddSpaceChargeField");

    if (pc.getCharge() == 0) {
        return;
    }

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!m_poisson_boundary_handler.bcs_set) {
        m_poisson_boundary_handler.definePhiBCs(Geom(0));
    }

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(n_rz_azimuthal_modes == 1,
                                     "Error: RZ electrostatic only implemented for a single mode");
#endif

    // Allocate fields for charge and potential
    const int num_levels = max_level + 1;
    Vector<std::unique_ptr<MultiFab> > rho(num_levels);
    Vector<std::unique_ptr<MultiFab> > rho_coarse(num_levels); // Used in order to interpolate between levels
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    // Use number of guard cells used for local deposition of rho
    const amrex::IntVect ng = guard_cells.ng_depos_rho;
    for (int lev = 0; lev <= max_level; lev++) {
        BoxArray nba = boxArray(lev);
        nba.surroundingNodes();
        rho[lev] = std::make_unique<MultiFab>(nba, DistributionMap(lev), 1, ng);
        rho[lev]->setVal(0.);
        phi[lev] = std::make_unique<MultiFab>(nba, DistributionMap(lev), 1, 1);
        phi[lev]->setVal(0.);
        if (lev > 0) {
            // For MR levels: allocated the coarsened version of rho
            BoxArray cba = nba;
            cba.coarsen(refRatio(lev-1));
            rho_coarse[lev] = std::make_unique<MultiFab>(cba, DistributionMap(lev), 1, ng);
            rho_coarse[lev]->setVal(0.);
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
    for (int lev = 0; lev <= max_level; lev++) {
        if (lev > 0) {
            if (charge_buf[lev]) {
                charge_buf[lev]->setVal(0.);
            }
        }
    }
    Warpx::SyncRho(rho, rho_coarse, charge_buf); // Apply filter, perform MPI exchange, interpolate across levels

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

/* Compute the potential `phi` by solving the Poisson equation with the
   simulation specific boundary conditions and boundary values, then add the
   E field due to that `phi` to `Efield_fp`.
*/
void
RelativisticExplicitES::AddBoundaryField (amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield, 
                                          amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield)
{
    WARPX_PROFILE("WarpX::AddBoundaryField");

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!m_poisson_boundary_handler.bcs_set) {
        m_poisson_boundary_handler.definePhiBCs(Geom(0));
    }

    // Allocate fields for charge and potential
    const int num_levels = max_level + 1;
    Vector<std::unique_ptr<MultiFab> > rho(num_levels);
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    // Use number of guard cells used for local deposition of rho
    const amrex::IntVect ng = guard_cells.ng_depos_rho;
    for (int lev = 0; lev <= max_level; lev++) {
        BoxArray nba = boxArray(lev);
        nba.surroundingNodes();
        rho[lev] = std::make_unique<MultiFab>(nba, DistributionMap(lev), 1, ng);
        rho[lev]->setVal(0.);
        phi[lev] = std::make_unique<MultiFab>(nba, DistributionMap(lev), 1, 1);
        phi[lev]->setVal(0.);
    }

    // Set the boundary potentials appropriately
    setPhiBC(phi);

    // beta is zero for boundaries
    const std::array<Real, 3> beta = {0._rt};

    // Compute the potential phi, by solving the Poisson equation
    computePhi( rho, phi, beta, self_fields_required_precision,
                self_fields_absolute_tolerance, self_fields_max_iters,
                self_fields_verbosity );

    // Compute the corresponding electric and magnetic field, from the potential phi.
    computeE( Efield_fp, phi, beta );
    computeB( Bfield_fp, phi, beta );
}

