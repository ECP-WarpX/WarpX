/* Copyright 2024 The WarpX Community
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

    // read input parameters
    const ParmParse pp_warpx("warpx");
    utils::parser::queryWithParser(pp_warpx, "semi_implicit_factor", m_C_SI);
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

    // The "sigma" multifabs stores the dressing of the Poisson equation. It
    // is a cell-centered multifab.
    WarpX::AllocInitMultiFab(
        sigma[lev], amrex::convert(ba, IntVect( AMREX_D_DECL(0,0,0) )),
        dm, ncomps, ngRho, lev, "sigma"
    );

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (ncomps == 1),
        "Implicit Darwin solver only support m = 0 azimuthal mode at present.");
#endif
}

void ImplicitDarwinSolver::ClearLevel (int lev)
{
    sigma[lev].reset();
    // for (int i = 0; i < 3; ++i) {
    //     sigma[lev][i].reset();
    // }
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
    sigma[lev]->setVal(1.0_rt);

    // GetChargeDensity returns a nodal multifab
    amrex::GpuArray<int, 3> nodal = {1, 1, 1};
    // sigma is a cell-centered array
    amrex::GpuArray<int, 3> const& cell_centered = {0, 0, 0};
    // The "coarsening is just 1 i.e. no coarsening"
    amrex::GpuArray<int, 3> const& coarsen = {1, 1, 1};

    // Below we set all the unused dimensions to have cell-centered values for
    // rho since these values will be interpolated onto a cell-centered grid
    // - if this is not done the Interp function returns nonsense values.
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_1D_Z)
    nodal[2] = 0;
#endif
#if defined(WARPX_DIM_1D_Z)
    nodal[1] = 0;
#endif

    auto& warpx = WarpX::GetInstance();
    auto& mypc = warpx.GetPartContainer();

    auto m_mult_factor = (
        m_C_SI * warpx.getdt(lev) * warpx.getdt(lev)
        / (16._rt * MathConst::pi * MathConst::pi * PhysConst::ep0)
    );

    // Loop over each species to calculate the Poisson equation dressing
    for (auto const& pc : mypc) {
        auto rho = pc->GetChargeDensity(lev, false);
        // Handle the parallel transfers of guard cells and
        // apply the filtering if requested - might not be needed...
        warpx.ApplyFilterandSumBoundaryRho(lev, lev, *rho, 0, rho->nComp());

        // multiply charge density by C_SI * dt**2 /4.0 q/(m \varepsilon_0)
        rho->mult(m_mult_factor * pc->getCharge() / pc->getMass());

        // update sigma
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*sigma[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            Array4<Real> const& sigma_arr = sigma[lev]->array(mfi);
            Array4<Real const> const& rho_arr = rho->const_array(mfi);

            // Loop over the cells and update the sigma field
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k){
                // Interpolate rho to cell-centered multifab to add to sigma.
                auto const rho_interp = ablastr::coarsen::sample::Interp(
                    rho_arr, nodal, cell_centered, coarsen, i, j, k, 0
                );
                sigma_arr(i, j, k, 0) += rho_interp;
            });

        }
    }
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
                "ImplicitDarwinSolver",
                "Max norm of rho is 0",
                ablastr::warn_manager::WarnPriority::low
        );
    }

    // auto& warpx = WarpX::GetInstance();
    // amrex::LPInfo info = LPInfo().setMaxCoarseningLevel(0);
    amrex::LPInfo info;

#if defined(AMREX_USE_EB) || defined(WARPX_DIM_RZ)
    // In the presence of EB or RZ: the solver assumes that the beam is
    // propagating along  one of the axes of the grid, i.e. that only *one*
    // of the components of `beta` is non-negligible.
    amrex::MLEBNodeLaplacian linop( {geom}, {rho.boxArray()}, {rho[lev].dm()}, info
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
    // linop.setScalars(0.0, 1.0);
    // linop.setBCoeffs(lev, amrex::GetArrOfConstPtrs(sigma[lev]) ); // for the non-axis-aligned solver
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
    // that can handle sigma in all directions
    amrex::MLNodeLaplacian linop(
        {geom}, {rho[lev]->boxArray()}, {rho[lev]->DistributionMap()}, info
        // {}, 1.0
    );
    linop.setSigma(lev, *sigma[lev]);
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

    // // create a vector to our fields, sorted by level
    // amrex::Vector<amrex::MultiFab*> sorted_rho;
    // amrex::Vector<amrex::MultiFab*> sorted_phi;
    // // for (int nlev = 0; nlev <= finest_level; ++nlev) {
    // sorted_rho.emplace_back(rho[lev].get());
    // sorted_phi.emplace_back(phi[lev].get());
    // // }

    // Solve Poisson equation at lev
    mlmg.solve(
        { phi[lev].get() }, { rho[lev].get() },   // GetVecOfPtrs(phi), GetVecOfConstPtrs(rho) //
        relative_tolerance, absolute_tolerance
    );

    // Run additional operations, such as calculation of the E field for embedded boundaries
    // if constexpr (!std::is_same<T_PostPhiCalculationFunctor, std::nullopt_t>::value)
    //     if (post_phi_calculation.has_value())
    //         post_phi_calculation.value()(mlmg, lev);

}
