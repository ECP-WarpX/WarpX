/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Remi Lehe, Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "FieldSolver/ElectrostaticSolvers/ElectrostaticSolver.H"
#include "Fluids/MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Python/callbacks.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <ablastr/fields/PoissonSolver.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MFIter.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>
#include <AMReX_MFInterp_C.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EBFabFactory.H>
#endif

#include <array>
#include <memory>
#include <string>

using namespace amrex;
using namespace warpx::fields;


void
WarpX::ComputeSpaceChargeField (bool const reset_fields)
{
    WARPX_PROFILE("WarpX::ComputeSpaceChargeField");
    if (reset_fields) {
        // Reset all E and B fields to 0, before calculating space-charge fields
        WARPX_PROFILE("WarpX::ComputeSpaceChargeField::reset_fields");
        for (int lev = 0; lev <= max_level; lev++) {
            for (int comp=0; comp<3; comp++) {
                Efield_fp[lev][comp]->setVal(0);
                Bfield_fp[lev][comp]->setVal(0);
            }
        }
    }

    m_electrostatic_solver->ComputeSpaceChargeField(Efield_fp, Bfield_fp);

    if (electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame ||
        electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic) {
        AddSpaceChargeFieldLabFrame();
    }
    else {
        // Loop over the species and add their space-charge contribution to E and B.
        // Note that the fields calculated here does not include the E field
        // due to simulation boundary potentials
        for (int ispecies=0; ispecies<mypc->nSpecies(); ispecies++){
            WarpXParticleContainer& species = mypc->GetParticleContainer(ispecies);
            if (species.initialize_self_fields ||
                (electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic)) {
                AddSpaceChargeField(species);
            }
        }

        // Add the field due to the boundary potentials
        if (m_boundary_potential_specified ||
                (electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic)){
            AddBoundaryField();
        }
    }
}

/* Compute the potential `phi` by solving the Poisson equation with the
   simulation specific boundary conditions and boundary values, then add the
   E field due to that `phi` to `Efield_fp`.
*/
void
WarpX::AddBoundaryField ()
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

void
WarpX::AddSpaceChargeField (WarpXParticleContainer& pc)
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
    SyncRho(rho, rho_coarse, charge_buf); // Apply filter, perform MPI exchange, interpolate across levels

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

void
WarpX::AddSpaceChargeFieldLabFrame ()
{
    WARPX_PROFILE("WarpX::AddSpaceChargeFieldLabFrame");

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!m_poisson_boundary_handler.bcs_set) {
        m_poisson_boundary_handler.definePhiBCs(Geom(0));
    }

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(n_rz_azimuthal_modes == 1,
                                     "Error: RZ electrostatic only implemented for a single mode");
#endif

    // Deposit particle charge density (source of Poisson solver)
    mypc->DepositCharge(rho_fp, 0.0_rt);
    if (do_fluid_species) {
        int const lev = 0;
        myfl->DepositCharge( lev, *rho_fp[lev] );
    }
    for (int lev = 0; lev <= max_level; lev++) {
        if (lev > 0) {
            if (charge_buf[lev]) {
                charge_buf[lev]->setVal(0.);
            }
        }
    }
    SyncRho(rho_fp, rho_cp, charge_buf); // Apply filter, perform MPI exchange, interpolate across levels
#ifndef WARPX_DIM_RZ
    for (int lev = 0; lev <= finestLevel(); lev++) {
        // Reflect density over PEC boundaries, if needed.
        ApplyRhofieldBoundary(lev, rho_fp[lev].get(), PatchType::fine);
    }
#endif

    // beta is zero in lab frame
    // Todo: use simpler finite difference form with beta=0
    const std::array<Real, 3> beta = {0._rt};

    // set the boundary potentials appropriately
    setPhiBC(phi_fp);

    // Compute the potential phi, by solving the Poisson equation
    if (IsPythonCallbackInstalled("poissonsolver")) {

        // Use the Python level solver (user specified)
        ExecutePythonCallback("poissonsolver");

    } else {

#if defined(WARPX_DIM_1D_Z)
        // Use the tridiag solver with 1D
        computePhiTriDiagonal(rho_fp, phi_fp);
#else
        // Use the AMREX MLMG or the FFT (IGF) solver otherwise
        computePhi(rho_fp, phi_fp, beta, self_fields_required_precision,
                   self_fields_absolute_tolerance, self_fields_max_iters,
                   self_fields_verbosity);
#endif

    }

    // Compute the electric field. Note that if an EB is used the electric
    // field will be calculated in the computePhi call.
    if (!m_eb_enabled) { computeE( Efield_fp, phi_fp, beta ); }
    else {
        if (IsPythonCallbackInstalled("poissonsolver")) { computeE(Efield_fp, phi_fp, beta); }
    }

    // Compute the magnetic field
    computeB( Bfield_fp, phi_fp, beta );
}

/* Compute the potential `phi` by solving the Poisson equation with `rho` as
   a source, assuming that the source moves at a constant speed \f$\vec{\beta}\f$.
   This uses the amrex solver.

   More specifically, this solves the equation
   \f[
       \vec{\nabla}^2 r \phi - (\vec{\beta}\cdot\vec{\nabla})^2 r \phi = -\frac{r \rho}{\epsilon_0}
   \f]

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
   \param[in] beta Represents the velocity of the source of `phi`
   \param[in] required_precision The relative convergence threshold for the MLMG solver
   \param[in] absolute_tolerance The absolute convergence threshold for the MLMG solver
   \param[in] max_iters The maximum number of iterations allowed for the MLMG solver
   \param[in] verbosity The verbosity setting for the MLMG solver
*/
void
WarpX::computePhi (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                   amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<Real, 3> const beta,
                   Real const required_precision,
                   Real absolute_tolerance,
                   int const max_iters,
                   int const verbosity) const {
    // create a vector to our fields, sorted by level
    amrex::Vector<amrex::MultiFab *> sorted_rho;
    amrex::Vector<amrex::MultiFab *> sorted_phi;
    for (int lev = 0; lev <= finest_level; ++lev) {
        sorted_rho.emplace_back(rho[lev].get());
        sorted_phi.emplace_back(phi[lev].get());
    }

    std::optional<ElectrostaticSolver::EBCalcEfromPhiPerLevel> post_phi_calculation;
#ifdef AMREX_USE_EB
    // TODO: double check no overhead occurs on "m_eb_enabled == false"
    std::optional<amrex::Vector<amrex::EBFArrayBoxFactory const *> > eb_farray_box_factory;
#else
    std::optional<amrex::Vector<amrex::FArrayBoxFactory const *> > const eb_farray_box_factory;
#endif
    if (m_eb_enabled)
    {
        // EB: use AMReX to directly calculate the electric field since with EB's the
        // simple finite difference scheme in WarpX::computeE sometimes fails
        if (electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame ||
            electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic)
        {
            // TODO: maybe make this a helper function or pass Efield_fp directly
            amrex::Vector<
                amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>
            > e_field;
            for (int lev = 0; lev <= finest_level; ++lev) {
                e_field.push_back(
#   if defined(WARPX_DIM_1D_Z)
                    amrex::Array<amrex::MultiFab*, 1>{
                        getFieldPointer(FieldType::Efield_fp, lev, 2)
                    }
#   elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    amrex::Array<amrex::MultiFab*, 2>{
                        getFieldPointer(FieldType::Efield_fp, lev, 0),
                        getFieldPointer(FieldType::Efield_fp, lev, 2)
                    }
#   elif defined(WARPX_DIM_3D)
                    amrex::Array<amrex::MultiFab *, 3>{
                        getFieldPointer(FieldType::Efield_fp, lev, 0),
                        getFieldPointer(FieldType::Efield_fp, lev, 1),
                        getFieldPointer(FieldType::Efield_fp, lev, 2)
                    }
#   endif
                );
            }
            post_phi_calculation = ElectrostaticSolver::EBCalcEfromPhiPerLevel(e_field);
        }

#ifdef AMREX_USE_EB
        amrex::Vector<
            amrex::EBFArrayBoxFactory const *
        > factories;
        for (int lev = 0; lev <= finest_level; ++lev) {
            factories.push_back(&WarpX::fieldEBFactory(lev));
        }
        eb_farray_box_factory = factories;
#endif
    }

    bool const is_solver_igf_on_lev0 =
        WarpX::poisson_solver_id == PoissonSolverAlgo::IntegratedGreenFunction;

    ablastr::fields::computePhi(
        sorted_rho,
        sorted_phi,
        beta,
        required_precision,
        absolute_tolerance,
        max_iters,
        verbosity,
        this->geom,
        this->dmap,
        this->grids,
        WarpX::grid_type,
        this->m_poisson_boundary_handler,
        is_solver_igf_on_lev0,
        m_eb_enabled,
        WarpX::do_single_precision_comms,
        this->ref_ratio,
        post_phi_calculation,
        gett_new(0),
        eb_farray_box_factory
    );

}

/* \brief Compute the potential by solving Poisson's equation with
          a 1D tridiagonal solve.

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
*/
void
WarpX::computePhiTriDiagonal (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi) const
{

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(max_level == 0,
        "The tridiagonal solver cannot be used with mesh refinement");

    const int lev = 0;

    const amrex::Real* dx = Geom(lev).CellSize();
    const amrex::Real xmin = Geom(lev).ProbLo(0);
    const amrex::Real xmax = Geom(lev).ProbHi(0);
    const int nx_full_domain = static_cast<int>( (xmax - xmin)/dx[0] + 0.5_rt );

    int nx_solve_min = 1;
    int nx_solve_max = nx_full_domain - 1;

    auto field_boundary_lo0 = WarpX::field_boundary_lo[0];
    auto field_boundary_hi0 = WarpX::field_boundary_hi[0];
    if (field_boundary_lo0 == FieldBoundaryType::Neumann || field_boundary_lo0 == FieldBoundaryType::Periodic) {
        // Neumann or periodic boundary condition
        // Solve for the point on the lower boundary
        nx_solve_min = 0;
    }
    if (field_boundary_hi0 == FieldBoundaryType::Neumann || field_boundary_hi0 == FieldBoundaryType::Periodic) {
        // Neumann or periodic boundary condition
        // Solve for the point on the upper boundary
        nx_solve_max = nx_full_domain;
    }

    // Create a 1-D MultiFab that covers all of x.
    // The tridiag solve will be done in this MultiFab and then copied out afterwards.
    const amrex::IntVect lo_full_domain(AMREX_D_DECL(0,0,0));
    const amrex::IntVect hi_full_domain(AMREX_D_DECL(nx_full_domain,0,0));
    const amrex::Box box_full_domain_node(lo_full_domain, hi_full_domain, amrex::IntVect::TheNodeVector());
    const BoxArray ba_full_domain_node(box_full_domain_node);
    const amrex::Vector<int> pmap = {0}; // The data will only be on processor 0
    const amrex::DistributionMapping dm_full_domain(pmap);

    // Put the data in the pinned arena since the tridiag solver will be done on the CPU, but have
    // the data readily accessible from the GPU.
    auto phi1d_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    auto zwork1d_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    auto rho1d_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));

    if (field_boundary_lo0 == FieldBoundaryType::PEC || field_boundary_hi0 == FieldBoundaryType::PEC) {
        // Copy from phi to get the boundary values
        phi1d_mf.ParallelCopy(*phi[lev], 0, 0, 1);
    }
    rho1d_mf.ParallelCopy(*rho[lev], 0, 0, 1);

    // Multiplier on the charge density
    const amrex::Real norm = dx[0]*dx[0]/PhysConst::ep0;
    rho1d_mf.mult(norm);

    // Use the MFIter loop since when parallel, only process zero has a FAB.
    // This skips the loop on all other processors.
    for (MFIter mfi(phi1d_mf); mfi.isValid(); ++mfi) {

        const auto& phi1d_arr = phi1d_mf[mfi].array();
        const auto& zwork1d_arr = zwork1d_mf[mfi].array();
        const auto& rho1d_arr = rho1d_mf[mfi].array();

        // The loops are always performed on the CPU

        amrex::Real diag = 2._rt;

        // The initial values depend on the boundary condition
        if (field_boundary_lo0 == FieldBoundaryType::PEC) {

            phi1d_arr(1,0,0) = (phi1d_arr(0,0,0) + rho1d_arr(1,0,0))/diag;

        } else if (field_boundary_lo0 == FieldBoundaryType::Neumann) {

            // Neumann boundary condition
            phi1d_arr(0,0,0) = rho1d_arr(0,0,0)/diag;

            zwork1d_arr(1,0,0) = 2._rt/diag;
            diag = 2._rt - zwork1d_arr(1,0,0);
            phi1d_arr(1,0,0) = (rho1d_arr(1,0,0) - (-1._rt)*phi1d_arr(1-1,0,0))/diag;

        } else if (field_boundary_lo0 == FieldBoundaryType::Periodic) {

            phi1d_arr(0,0,0) = rho1d_arr(0,0,0)/diag;

            zwork1d_arr(1,0,0) = 1._rt/diag;
            diag = 2._rt - zwork1d_arr(1,0,0);
            phi1d_arr(1,0,0) = (rho1d_arr(1,0,0) - (-1._rt)*phi1d_arr(1-1,0,0))/diag;

        }

        // Loop upward, calculating the Gaussian elimination multipliers and right hand sides
        for (int i_up = 2 ; i_up < nx_solve_max ; i_up++) {

            zwork1d_arr(i_up,0,0) = 1._rt/diag;
            diag = 2._rt - zwork1d_arr(i_up,0,0);
            phi1d_arr(i_up,0,0) = (rho1d_arr(i_up,0,0) - (-1._rt)*phi1d_arr(i_up-1,0,0))/diag;

        }

        // The last value depend on the boundary condition
        amrex::Real zwork_product = 1.; // Needed for parallel boundaries
        if (field_boundary_hi0 == FieldBoundaryType::PEC) {

            int const nxm1 = nx_full_domain - 1;
            zwork1d_arr(nxm1,0,0) = 1._rt/diag;
            diag = 2._rt - zwork1d_arr(nxm1,0,0);
            phi1d_arr(nxm1,0,0) = (phi1d_arr(nxm1+1,0,0) + rho1d_arr(nxm1,0,0) - (-1._rt)*phi1d_arr(nxm1-1,0,0))/diag;

        } else if (field_boundary_hi0 == FieldBoundaryType::Neumann) {

            // Neumann boundary condition
            zwork1d_arr(nx_full_domain,0,0) = 1._rt/diag;
            diag = 2._rt - 2._rt*zwork1d_arr(nx_full_domain,0,0);
            if (diag == 0._rt) {
                // This happens if the lower boundary is also Neumann.
                // It this case, the potential is relative to an arbitrary constant,
                // so set the upper boundary to zero to force a value.
                phi1d_arr(nx_full_domain,0,0) = 0.;
            } else {
                phi1d_arr(nx_full_domain,0,0) = (rho1d_arr(nx_full_domain,0,0) - (-1._rt)*phi1d_arr(nx_full_domain-1,0,0))/diag;
            }

        } else if (field_boundary_hi0 == FieldBoundaryType::Periodic) {

            zwork1d_arr(nx_full_domain,0,0) = 1._rt/diag;

            for (int i = 1 ; i <= nx_full_domain ; i++) {
                zwork_product *= zwork1d_arr(i,0,0);
            }

            diag = 2._rt - zwork1d_arr(nx_full_domain,0,0) - zwork_product;
            // Note that rho1d_arr(0,0,0) is used to ensure that the same value is used
            // on both boundaries.
            phi1d_arr(nx_full_domain,0,0) = (rho1d_arr(0,0,0) - (-1._rt)*phi1d_arr(nx_full_domain-1,0,0))/diag;

        }

        // Loop downward to calculate the phi
        if (field_boundary_lo0 == FieldBoundaryType::Periodic) {

            // With periodic, the right hand column adds an extra term for all rows
            for (int i_down = nx_full_domain-1 ; i_down >= 0 ; i_down--) {
                zwork_product /= zwork1d_arr(i_down+1,0,0);
                phi1d_arr(i_down,0,0) = phi1d_arr(i_down,0,0) + zwork1d_arr(i_down+1,0,0)*phi1d_arr(i_down+1,0,0) + zwork_product*phi1d_arr(nx_full_domain,0,0);
            }

        } else {

            for (int i_down = nx_solve_max-1 ; i_down >= nx_solve_min ; i_down--) {
                phi1d_arr(i_down,0,0) = phi1d_arr(i_down,0,0) + zwork1d_arr(i_down+1,0,0)*phi1d_arr(i_down+1,0,0);
            }

        }

    }

    // Copy phi1d to phi
    phi[lev]->ParallelCopy(phi1d_mf, 0, 0, 1);
}
