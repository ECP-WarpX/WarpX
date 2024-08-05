/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "SemiImplicitSolver.H"

using namespace amrex;

SemiImplicitSolver::SemiImplicitSolver ( int nlevs_max )
{
    AllocateMFs(nlevs_max);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (nlevs_max == 1),
        "Semi-implicit electrostatic solver only supports one level at present"
    );

    // read input parameters
    const ParmParse pp_warpx("warpx");
    utils::parser::queryWithParser(pp_warpx, "semi_implicit_factor", m_C_SI);
}

void SemiImplicitSolver::AllocateMFs (int nlevs_max)
{
    sigma.resize(nlevs_max);
}

void SemiImplicitSolver::AllocateLevelMFs (
    int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
    int ncomps, const amrex::IntVect& ngRho, const amrex::IntVect& rho_nodal_flag)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (rho_nodal_flag == IntVect::TheNodeVector()),
        "Semi-implicit electrostatic solver only supports a nodal charge density field"
    );

    // The "sigma" multifab stores the dressing of the Poisson equation. It
    // is a cell-centered multifab.
    WarpX::AllocInitMultiFab(
        sigma[lev], amrex::convert(ba, IntVect( AMREX_D_DECL(0,0,0) )),
        dm, ncomps, ngRho, lev, "sigma"
    );

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (ncomps == 1),
        "Semi-implicit electrostatic solver only supports m = 0 azimuthal mode at present.");
#endif
}

void SemiImplicitSolver::ClearLevel (int lev)
{
    sigma[lev].reset();
}

void SemiImplicitSolver::ComputeSigma () {

    int const lev = 0;
    sigma[lev]->setVal(1.0_rt);

    // sigma is a cell-centered array
    amrex::GpuArray<int, 3> const cell_centered = {0, 0, 0};
    // The "coarsening is just 1 i.e. no coarsening"
    amrex::GpuArray<int, 3> const coarsen = {1, 1, 1};

    // GetChargeDensity returns a nodal multifab
    // Below we set all the unused dimensions to have cell-centered values for
    // rho since these values will be interpolated onto a cell-centered grid
    // - if this is not done the Interp function returns nonsense values.
#if defined(WARPX_DIM_3D)
    amrex::GpuArray<int, 3> const nodal = {1, 1, 1};
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    amrex::GpuArray<int, 3> const nodal = {1, 1, 0};
#elif defined(WARPX_DIM_1D_Z)
    amrex::GpuArray<int, 3> const nodal = {1, 0, 0};
#endif

    auto& warpx = WarpX::GetInstance();
    auto& mypc = warpx.GetPartContainer();

    // The semi-implicit dielectric function is given by
    // \varepsilon_{SI} = \varepsilon * (1 + \sum_{i in species} C_{SI}*(w_pi * dt)^2/4)
    // Note the use of the plasma frequency in rad/s (not Hz) and the factor of 1/4,
    // these choices make it so that C_SI = 1 is the marginal stability threshold.
    auto mult_factor = (
        m_C_SI * warpx.getdt(lev) * warpx.getdt(lev) / (4._rt * PhysConst::ep0)
    );

    // Loop over each species to calculate the Poisson equation dressing
    for (auto const& pc : mypc) {
        // grab the charge density for this species
        auto rho = pc->GetChargeDensity(lev, false);

        // Handle the parallel transfer of guard cells and apply filtering
        warpx.ApplyFilterandSumBoundaryRho(lev, lev, *rho, 0, rho->nComp());

        // get multiplication factor for this species
        auto const mult_factor_pc = mult_factor * pc->getCharge() / pc->getMass();

        // update sigma
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*sigma[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            Array4<Real> const& sigma_arr = sigma[lev]->array(mfi);
            Array4<Real const> const& rho_arr = rho->const_array(mfi);

            // Loop over the cells and update the sigma field
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k){
                // Interpolate rho to cell-centered value
                auto const rho_cc = ablastr::coarsen::sample::Interp(
                    rho_arr, nodal, cell_centered, coarsen, i, j, k, 0
                );
                // add species term to sigma:
                // C_SI * w_p^2 * dt^2 / 4 = C_SI / 4 * q*rho/(m*eps0) * dt^2
                sigma_arr(i, j, k, 0) += mult_factor_pc * rho_cc;
            });

        }
    }
}
