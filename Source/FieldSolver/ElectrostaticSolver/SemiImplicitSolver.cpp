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
