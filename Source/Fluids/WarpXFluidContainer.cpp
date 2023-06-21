/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote
 * Jean-Luc Vay, Luca Fedeli, Maxence Thevenet
 * Michael Rowan, Remi Lehe, Revathi Jambunathan
 * Weiqun Zhang, Yinjian Zhao, levinem
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ablastr/coarsen/sample.H"

#include "WarpXFluidContainer.H"
#include "WarpX.H"
#include <ablastr/utils/Communication.H>
using namespace ablastr::utils::communication;
using namespace amrex;

WarpXFluidContainer::WarpXFluidContainer(int nlevs_max, int ispecies, const std::string &name)
{
    species_id = ispecies;
    species_name = name;

    plasma_injector = std::make_unique<PlasmaInjector>(species_id, species_name);
    physical_species = plasma_injector->getPhysicalSpecies();
    charge = plasma_injector->getCharge();
    mass = plasma_injector->getMass();

    ReadParameters();

    // Resize the list of MultiFabs for the right number of levels
    N.resize(nlevs_max);
    NU.resize(nlevs_max);
}

void WarpXFluidContainer::ReadParameters()
{
    static bool initialized = false;
    if (!initialized)
    {
        const ParmParse pp_species_name(species_name);
        pp_species_name.query("do_not_deposit", do_not_deposit);
        pp_species_name.query("do_not_gather", do_not_gather);
        pp_species_name.query("do_not_push", do_not_push);
        initialized = true;
    }
}

void WarpXFluidContainer::AllocateLevelMFs(int lev, const BoxArray &ba, const DistributionMapping &dm)
{
    int ncomps = 1;
    amrex::IntVect nguards = {AMREX_D_DECL(2, 2, 2)};

    auto &warpx = WarpX::GetInstance();

    // set human-readable tag for each MultiFab
    auto const tag = [lev](std::string tagname)
    {
        tagname.append("[l=").append(std::to_string(lev)).append("]");
        return tagname;
    };

    warpx.AllocInitMultiFab(N[lev], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, tag("fluid density"), 0.0_rt);

    warpx.AllocInitMultiFab(NU[lev][0], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, tag("fluid momentum density [x]"), 0.0_rt);
    warpx.AllocInitMultiFab(NU[lev][1], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, tag("fluid momentum density [y]"), 0.0_rt);
    warpx.AllocInitMultiFab(NU[lev][2], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, tag("fluid momentum density [z]"), 0.0_rt);
}

void WarpXFluidContainer::InitData(int lev)
{
    // Extract objects that give the initial density and momentum
    InjectorDensity *inj_rho = plasma_injector->getInjectorDensity();
    InjectorMomentum *inj_mom = plasma_injector->getInjectorMomentum();

    // Extract grid geometry properties
    WarpX &warpx = WarpX::GetInstance();
    const amrex::Geometry &geom = warpx.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    // Loop through cells and initialize their value
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        amrex::Box const &tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());
        amrex::Array4<Real> const &N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> const &NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> const &NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> const &NUz_arr = NU[lev][2]->array(mfi);

        amrex::ParallelFor(tile_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                           {
#if defined(WARPX_DIM_3D)
                               amrex::Real x = problo[0] + i * dx[0];
                               amrex::Real y = problo[1] + j * dx[1];
                               amrex::Real z = problo[2] + k * dx[2];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                amrex::Real x = problo[0] + i * dx[0];
                amrex::Real y = 0.0_rt;
                amrex::Real z = problo[1] + j * dx[1];
#else
                amrex::Real x = 0.0_rt;
                amrex::Real y = 0.0_rt;
                amrex::Real z = problo[0] + i * dx[0];
#endif

                               amrex::Real n = inj_rho->getDensity(x, y, z);
                               auto u = inj_mom->getBulkMomentum(x, y, z);

                               N_arr(i, j, k) = n;
                               NUx_arr(i, j, k) = n * u.x;
                               NUy_arr(i, j, k) = n * u.y;
                               NUz_arr(i, j, k) = n * u.z; });
    }

    // Fill guard cells
    const amrex::Periodicity &period = geom.periodicity();
    FillBoundary(*N[lev], N[lev]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*N[lev], NU[0][lev]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*N[lev], NU[1][lev]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*N[lev], NU[2][lev]->nGrowVect(), WarpX::do_single_precision_comms, period);
}

void WarpXFluidContainer::Evolve(int lev,
                                 const amrex::MultiFab &Ex, const amrex::MultiFab &Ey, const amrex::MultiFab &Ez,
                                 const amrex::MultiFab &Bx, const amrex::MultiFab &By, const amrex::MultiFab &Bz,
                                 amrex::MultiFab &jx, amrex::MultiFab &jy, amrex::MultiFab &jz,
                                 amrex::MultiFab *cjx, amrex::MultiFab *cjy, amrex::MultiFab *cjz,
                                 amrex::MultiFab *rho, amrex::MultiFab *crho,
                                 const amrex::MultiFab *cEx, const amrex::MultiFab *cEy, const amrex::MultiFab *cEz,
                                 const amrex::MultiFab *cBx, const amrex::MultiFab *cBy, const amrex::MultiFab *cBz,
                                 amrex::Real t, amrex::Real dt, DtType a_dt_type, bool skip_deposition)
{
    // Gather E&B fields to each node
    // TODO

    // Update N and NU at each node
    // TODO

    // Deposit J to the simulation mesh
    if (!skip_deposition)
    {
        DepositCurrent(lev, jx, jy, jz);
    }
}

void WarpXFluidContainer::DepositCurrent(int lev,
                                         amrex::MultiFab &jx, amrex::MultiFab &jy, amrex::MultiFab &jz)
{

    const amrex::ParticleReal inv_clight_sq = 1.0_prt / PhysConst::c / PhysConst::c;
    charge = plasma_injector->getCharge();

// Loop over grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        amrex::Box const &tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());
        amrex::Array4<Real> const &N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> const &NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> const &NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> const &NUz_arr = NU[lev][2]->array(mfi);

        amrex::Array4<amrex::Real> jx_arr = jx.array(mfi);
        amrex::Array4<amrex::Real> jy_arr = jy.array(mfi);
        amrex::Array4<amrex::Real> jz_arr = jz.array(mfi);

        // Loop over cells (ParallelFor)
        amrex::ParallelFor(tile_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                           {
                               // Interpolate N/NU from nodes to the simulation mesh (typically Yee mesh)
                               amrex::GpuArray<int, 3U> sf = {AMREX_D_DECL(1, 1, 1)};
                               auto N_CC = ablastr::coarsen::sample::Interp(N_arr, sf, sf, sf, i, j, k, 1);
                               auto NUx_CC = ablastr::coarsen::sample::Interp(NUx_arr, sf, sf, sf, i, j, k, 1);
                               auto NUy_CC = ablastr::coarsen::sample::Interp(NUy_arr, sf, sf, sf, i, j, k, 1);
                               auto NUz_CC = ablastr::coarsen::sample::Interp(NUz_arr, sf, sf, sf, i, j, k, 1);

                               // Calculate J from fluid and add it to jx/jy/jz
                               auto gamma = std::sqrt(N_CC * N_CC + (NUx_CC * NUx_CC + NUy_CC * NUy_CC + NUz_CC * NUz_CC)*inv_clight_sq) / N_CC;
                               jx_arr(i, j, k) += charge * (NUx_CC / gamma);
                               jy_arr(i, j, k) += charge * (NUy_CC / gamma);
                               jz_arr(i, j, k) += charge * (NUz_CC / gamma); });
    }
}
