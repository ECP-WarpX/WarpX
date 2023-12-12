/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ImplicitDarwinSolver.H"

using namespace amrex;

ImplicitDarwinSolver::ImplicitDarwinSolver ( int nlevs_max )
{
    AllocateMFs(nlevs_max);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (nlevs_max == 1),
        "Implicit Darwin solver only support one level at present"
    );
}

void ImplicitDarwinSolver::AllocateMFs (int nlevs_max)
{
    sigma.resize(nlevs_max);
}

void ImplicitDarwinSolver::AllocateLevelMFs (
    int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
    int ncomps, const amrex::IntVect& ngRho, const amrex::IntVect& rho_nodal_flag)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (rho_nodal_flag == IntVect::TheNodeVector()),
        "Implicit Darwin solver only support a nodal charge density field at present"
    );

    // sigma[lev][0] = MultiFab(
    //     amrex::convert(ba, IntVect( AMREX_D_DECL(0,0,0) )),
    //     dm, ncomps, ngRho
    // );
    // sigma[lev][0].setVal(0.0_rt);
    // sigma[lev][1] = MultiFab(
    //     amrex::convert(ba, IntVect( AMREX_D_DECL(0,0,0) )),
    //     dm, ncomps, ngRho
    // );
    // sigma[lev][1].setVal(0.0_rt);

    // The "sigma" multifabs stores the dressing of the Poisson equation that
    // eliminates plasma frequency effects.
    WarpX::AllocInitMultiFab(
        sigma[lev][0], amrex::convert(ba, IntVect( AMREX_D_DECL(0,0,0) )),
        dm, ncomps, ngRho, lev, "sigma[x]", 1.0_rt
    );
    WarpX::AllocInitMultiFab(
        sigma[lev][1], amrex::convert(ba, IntVect( AMREX_D_DECL(0,0,0) )),
        dm, ncomps, ngRho, lev, "sigma[y]", 1.0_rt
    );
    // WarpX::AllocInitMultiFab(
    //     sigma[lev][2], amrex::convert(ba, IntVect( AMREX_D_DECL(0,0,0) )),
    //     dm, ncomps, ngRho, lev, "sigma[z]", 0.0_rt
    // );
    // sigma[lev][1] = std::make_unique<amrex::MultiFab>(
    //     sigma[lev][0], amrex::make_alias, 0, ncomps
    // );
    // sigma[lev][2] = std::make_unique<amrex::MultiFab>(
    //     sigma[lev][0], amrex::make_alias, 0, ncomps
    // );

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (ncomps == 1),
        "Implicit Darwin solver only support m = 0 azimuthal mode at present.");
#endif
}

void ImplicitDarwinSolver::ClearLevel (int lev)
{
    // sigma[lev].reset();
    for (int i = 0; i < 3; ++i) {
        sigma[lev][i].reset();
    }
}

void
ImplicitDarwinSolver::AddSpaceChargeField (
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>& rhofield,
    amrex::Vector<std::unique_ptr<amrex::MultiFab>>& phifield,
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
    Real const required_precision, Real absolute_tolerance,
    int const max_iters, int const verbosity
)
{
    // Reset all E fields to 0, before calculating space-charge fields
    WARPX_PROFILE("WarpX::ComputeSpaceChargeField::reset_fields");
    for (int lev = 0; lev <= 0; lev++) {
        for (int comp=0; comp<3; comp++) {
            Efield[lev][comp]->setVal(0);
            // Bfield_fp[lev][comp]->setVal(0);
        }
    }

    WARPX_PROFILE("ImplicitDarwinSolver::AddSpaceChargeField");
    auto& warpx = WarpX::GetInstance();

    // Store the boundary conditions for the field solver if they haven't been
    // stored yet
    if (!warpx.m_poisson_boundary_handler.bcs_set)
        warpx.m_poisson_boundary_handler.definePhiBCs(warpx.Geom(0));

    // Deposit particle charge density (source of Poisson solver)
    auto& mypc = warpx.GetPartContainer();
    mypc.DepositCharge(rhofield, 0.0_rt);
    // if (warpx.do_fluid_species) {
    //     int const lev = 0;
    //     myfl->DepositCharge( lev, *rho_fp[lev] );
    // }

    // Reflect density over PEC boundaries, if needed.
    // filter (if used), exchange guard cells, interpolate across MR levels
    // and apply boundary conditions
    warpx.SyncCurrentAndRho();

    // sigma dresses the Poisson equation to filter out plasma frequency effects
    ComputeSigma();

#if defined(AMREX_USE_EB)
    std::optional<amrex::Vector<amrex::EBFArrayBoxFactory const *> > eb_farray_box_factory;
    amrex::Vector<
        amrex::EBFArrayBoxFactory const *
    > factories;
    for (int lev = 0; lev <= finest_level; ++lev) {
        factories.push_back(&warpx.fieldEBFactory(lev));
    }
    eb_farray_box_factory = factories;
#else
    const std::optional<ElectrostaticSolver::EBCalcEfromPhiPerLevel> post_phi_calculation;
    const std::optional<amrex::Vector<amrex::FArrayBoxFactory const *> > eb_farray_box_factory;
#endif

    // Use the AMREX MLMG solver
    ComputePhi(
        rhofield, phifield, required_precision,
        absolute_tolerance, max_iters, verbosity,
        warpx.Geom(0),
        warpx.m_poisson_boundary_handler,
        eb_farray_box_factory
    );

    // Compute the electric field. Note that if an EB is used the electric
    // field will be calculated in the computePhi call above.
#ifndef AMREX_USE_EB
    // beta is the moving frame velocity i.e. zero in lab frame
    const std::array<Real, 3> beta = {0._rt};

    warpx.computeE( Efield, phifield, beta );
#endif

//     // Compute the magnetic field
//     computeB( Bfield_fp, phi_fp, beta );
}

void ImplicitDarwinSolver::ComputeSigma () {

    int const lev = 0;

    auto& warpx = WarpX::GetInstance();
    auto& mypc = warpx.GetPartContainer();
    auto& pc_electron = mypc.GetParticleContainerFromName("electrons");

    // GetChargeDensity returns a cell-centered multifab
    auto rho_electron = pc_electron.GetChargeDensity(lev, false);

    // multiply charge density by C_SI * dt**2 /4.0 q/(m \varepsilon_0)
    rho_electron->mult(
        (2.0_rt / 4.0_rt) * warpx.getdt(lev)*warpx.getdt(lev)
        * pc_electron.getCharge()
        / (pc_electron.getMass())
    );

    // set sigma multifabs to 1 then add rho_electron to it
    for (int ii = 0; ii < AMREX_SPACEDIM; ii++)
    {
        sigma[lev][ii]->setVal(1.0_rt);
        MultiFab::Add(
            *sigma[lev][ii], *rho_electron, 0, 0, 1,
            sigma[lev][ii]->nGrowVect()
        );
    }


    // // Temporary cell-centered, single-component MultiFab for storing particles per cell.
    // amrex::MultiFab ppc_mf(warpx.boxArray(lev), warpx.DistributionMap(lev), 1, 1);
    // // Set value to 0, and increment the value in each cell with ppc.
    // ppc_mf.setVal(0._rt);
    // // Compute ppc which includes a summation over all species.
    // pc_electron.Increment(ppc_mf, lev);
}

template<typename T_BoundaryHandler, typename T_FArrayBoxFactory>
void ImplicitDarwinSolver::ComputePhi (
    const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
    amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
    Real const relative_tolerance,
    Real absolute_tolerance,
    int const max_iters,
    int const verbosity,
    amrex::Geometry const geom,
    T_BoundaryHandler const boundary_handler,
    std::optional<amrex::Vector<T_FArrayBoxFactory const *> > eb_farray_box_factory ) const
{
    using namespace amrex::literals;

    int const lev = 0;

    // scale rho appropriately; also determine if rho is zero everywhere
    amrex::Real max_norm_b = 0.0;
    rho[lev]->mult(-1._rt/ablastr::constant::SI::ep0);  // TODO: when do we "un-multiply" this? We need to document this side-effect!
    max_norm_b = amrex::max(max_norm_b, rho[lev]->norm0());

    amrex::ParallelDescriptor::ReduceRealMax(max_norm_b);

    const bool always_use_bnorm = (max_norm_b > 0);
    if (!always_use_bnorm) {
        if (absolute_tolerance == 0.0) absolute_tolerance = amrex::Real(1e-6);
        ablastr::warn_manager::WMRecordWarning(
                "ElectrostaticSolver",
                "Max norm of rho is 0",
                ablastr::warn_manager::WarnPriority::low
        );
    }

    amrex::LPInfo info;

#if defined(AMREX_USE_EB) || defined(WARPX_DIM_RZ)
    // In the presence of EB or RZ: the solver assumes that the beam is
    // propagating along  one of the axes of the grid, i.e. that only *one*
    // of the components of `beta` is non-negligible.
    amrex::MLABecLaplacian linop( {geom}, {rho.boxArray()}, {rho[lev].dm()}, info
#if defined(AMREX_USE_EB)
        , {eb_farray_box_factory.value()[lev]}
#endif
    );

    // Note: this assumes that the beam is propagating along
    // one of the axes of the grid, i.e. that only *one* of the
    // components of `beta` is non-negligible. // we use this
#if defined(WARPX_DIM_RZ)
    linop.setSigma({0._rt, 1._rt-beta_solver[1]*beta_solver[1]});
#else
    linop.setSigma({AMREX_D_DECL(
        1._rt-beta_solver[0]*beta_solver[0],
        1._rt-beta_solver[1]*beta_solver[1],
        1._rt-beta_solver[2]*beta_solver[2])});
#endif

#if defined(AMREX_USE_EB)
    // if the EB potential only depends on time, the potential can be passed
    // as a float instead of a callable
    if (boundary_handler.phi_EB_only_t) {
        linop.setEBDirichlet(boundary_handler.potential_eb_t(current_time.value()));
    }
    else
        linop.setEBDirichlet(boundary_handler.getPhiEB(current_time.value()));
#endif
#else
    // In the absence of EB and RZ: use a more generic solver
    // that can handle beams propagating in any direction
    amrex::MLABecLaplacian linop(
        {geom}, {rho[lev]->boxArray()}, {rho[lev]->DistributionMap()}, info
    );
    linop.setScalars(0.0, 1.0);
    linop.setBCoeffs(lev, amrex::GetArrOfConstPtrs(sigma[lev]) ); // for the non-axis-aligned solver
#endif

    // Solve the Poisson equation
    linop.setDomainBC( boundary_handler.lobc, boundary_handler.hibc );
#ifdef WARPX_DIM_RZ
    linop.setRZ(true);
#endif
    amrex::MLMG mlmg(linop); // actual solver defined here
    mlmg.setVerbose(verbosity);
    mlmg.setMaxIter(max_iters);
    mlmg.setAlwaysUseBNorm(always_use_bnorm);

    // Solve Poisson equation at lev
    mlmg.solve(
        GetVecOfPtrs(phi), GetVecOfConstPtrs(rho),
        relative_tolerance, absolute_tolerance
    );

    // Run additional operations, such as calculation of the E field for embedded boundaries
    // if constexpr (!std::is_same<T_PostPhiCalculationFunctor, std::nullopt_t>::value)
    //     if (post_phi_calculation.has_value())
    //         post_phi_calculation.value()(mlmg, lev);

}
