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
#include "Particles/Pusher/UpdateMomentumHigueraCary.H"

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

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
                NUz_arr(i, j, k) = n * u.z;
            }
        );
    }

    // Fill guard cells
    const amrex::Periodicity &period = geom.periodicity();
    FillBoundary(*N[lev], N[lev]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][0], NU[lev][0]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][1], NU[lev][1]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][2], NU[lev][2]->nGrowVect(), WarpX::do_single_precision_comms, period);
}


void WarpXFluidContainer::Evolve(
    int lev,
    const amrex::MultiFab &Ex, const amrex::MultiFab &Ey, const amrex::MultiFab &Ez,
    const amrex::MultiFab &Bx, const amrex::MultiFab &By, const amrex::MultiFab &Bz,
    amrex::MultiFab &jx, amrex::MultiFab &jy, amrex::MultiFab &jz,
    amrex::MultiFab *cjx, amrex::MultiFab *cjy, amrex::MultiFab *cjz,
    amrex::MultiFab *rho, amrex::MultiFab *crho,
    const amrex::MultiFab *cEx, const amrex::MultiFab *cEy, const amrex::MultiFab *cEz,
    const amrex::MultiFab *cBx, const amrex::MultiFab *cBy, const amrex::MultiFab *cBz,
    amrex::Real t, amrex::Real dt, DtType a_dt_type, bool skip_deposition)
{
    GatherAndPush(lev, Ex, Ey, Ez, Bx, By, Bz);

    // Deposit J to the simulation mesh
    if (!skip_deposition)
    {
        DepositCurrent(lev, jx, jy, jz);
    }
}


void WarpXFluidContainer::GatherAndPush (
    int lev,
    const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
    const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz)
{

    WarpX &warpx = WarpX::GetInstance();
    const amrex::Real q = getCharge();
    const amrex::Real m = getMass();
    const Real dt = warpx.getdt(lev);

    // Temporary nodal Fields
    amrex::MultiFab tmp_Ex_Nodal(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_Ey_Nodal(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_Ez_Nodal(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_Bx_Nodal(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_By_Nodal(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_Bz_Nodal(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);

   // Prepare interpolation of current components to cell center
    // The arrays below store the index type (staggering) of each MultiFab, with the third
    // component set to zero in the two-dimensional case.
    auto Nodal_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Ex_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Ey_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Ez_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Bx_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto By_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Bz_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        Nodal_type[i] = N[lev]->ixType()[i];
        Ex_CC_type[i] = Ex.ixType()[i];
        Ey_CC_type[i] = Ey.ixType()[i];
        Ez_CC_type[i] = Ez.ixType()[i];
        Bx_CC_type[i] = Bx.ixType()[i];
        By_CC_type[i] = By.ixType()[i];
        Bz_CC_type[i] = Bz.ixType()[i];
    }

    // MFIter loop
    //    ParallelFor loop (over nodes)
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        amrex::Box const &tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());

        amrex::Array4<amrex::Real> tmp_Ex_Nodal_arr = tmp_Ex_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Ey_Nodal_arr = tmp_Ey_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Ez_Nodal_arr = tmp_Ez_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Bx_Nodal_arr = tmp_Bx_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_By_Nodal_arr = tmp_By_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Bz_Nodal_arr = tmp_Bz_Nodal.array(mfi);

        amrex::Array4<const amrex::Real> const& Ex_arr = Ex.array(mfi);
        amrex::Array4<const amrex::Real> const& Ey_arr = Ey.array(mfi);
        amrex::Array4<const amrex::Real> const& Ez_arr = Ez.array(mfi);
        amrex::Array4<const amrex::Real> const& Bx_arr = Bx.array(mfi);
        amrex::Array4<const amrex::Real> const& By_arr = By.array(mfi);
        amrex::Array4<const amrex::Real> const& Bz_arr = Bz.array(mfi);

        // Here, we do not perform any coarsening.
        amrex::GpuArray<int, 3U> coarsening_ratio = {AMREX_D_DECL(1, 1, 1)};

  // Loop over cells (ParallelFor)
        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // Interpolate fields from CC to Nodal points
                amrex::Real Ex_Nodal = ablastr::coarsen::sample::Interp(Ex_arr,
                    Ex_CC_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Ey_Nodal = ablastr::coarsen::sample::Interp(Ey_arr,
                    Ey_CC_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Ez_Nodal = ablastr::coarsen::sample::Interp(Ez_arr,
                    Ez_CC_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Bx_Nodal = ablastr::coarsen::sample::Interp(Bx_arr,
                    Bx_CC_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real By_Nodal = ablastr::coarsen::sample::Interp(By_arr,
                    By_CC_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Bz_Nodal = ablastr::coarsen::sample::Interp(Bz_arr,
                    Bz_CC_type, Nodal_type, coarsening_ratio, i, j, k, 0);

                tmp_Ex_Nodal_arr(i, j, k) += Ex_Nodal;
                tmp_Ey_Nodal_arr(i, j, k) += Ey_Nodal;  
                tmp_Ez_Nodal_arr(i, j, k) += Ez_Nodal;
                tmp_Bx_Nodal_arr(i, j, k) += Bx_Nodal;
                tmp_By_Nodal_arr(i, j, k) += By_Nodal;
                tmp_Bz_Nodal_arr(i, j, k) += Bz_Nodal;

            }
        );
    }


    // MFIter loop
    //    ParallelFor loop (over nodes)
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        amrex::Box const &tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());

        amrex::Array4<Real> const &N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> NUz_arr = NU[lev][2]->array(mfi);

        amrex::Array4<amrex::Real> tmp_Ex_Nodal_arr = tmp_Ex_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Ey_Nodal_arr = tmp_Ey_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Ez_Nodal_arr = tmp_Ez_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Bx_Nodal_arr = tmp_Bx_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_By_Nodal_arr = tmp_By_Nodal.array(mfi);
        amrex::Array4<amrex::Real> tmp_Bz_Nodal_arr = tmp_Bz_Nodal.array(mfi);

        // Here, we do not perform any coarsening.
        amrex::GpuArray<int, 3U> coarsening_ratio = {AMREX_D_DECL(1, 1, 1)};

        // Loop over cells (ParallelFor)
        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // Calculate J from fluid and add it to jx/jy/jz
                auto tmp_Ux = (NUx_arr(i, j, k) / N_arr(i,j,k));
                auto tmp_Uy = (NUy_arr(i, j, k) / N_arr(i,j,k));
                auto tmp_Uz = (NUz_arr(i, j, k) / N_arr(i,j,k));

                // Push the fluid elements momentum 
                UpdateMomentumHigueraCary(tmp_Ux, tmp_Uy, tmp_Uz,
                    tmp_Ex_Nodal_arr(i, j, k), tmp_Ey_Nodal_arr(i, j, k), tmp_Ez_Nodal_arr(i, j, k), 
                    tmp_Bx_Nodal_arr(i, j, k), tmp_By_Nodal_arr(i, j, k), tmp_Bz_Nodal_arr(i, j, k), q, m, dt );

                // Calculate NU
                NUx_arr(i,j,k) = N_arr(i,j,k)*tmp_Ux;
                NUy_arr(i,j,k) = N_arr(i,j,k)*tmp_Uy;
                NUz_arr(i,j,k) = N_arr(i,j,k)*tmp_Uz;
                
            }
        );
    }




    //        - Interpolate E and B to the nodal grid
    //          (store the values in local variables of type `amrex::Real`,
    //           similar to `jy_CC` in `DepositCurrent)
    //        - Use these values to update `U`
    //          by calling UpdateMomentumHigueraCary

}

void WarpXFluidContainer::DepositCurrent(
    int lev,
    amrex::MultiFab &jx, amrex::MultiFab &jy, amrex::MultiFab &jz)
{

    // Temporary nodal currents
    amrex::MultiFab tmp_jx_fluid(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_jy_fluid(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_jz_fluid(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);

    const amrex::Real inv_clight_sq = 1.0_prt / PhysConst::c / PhysConst::c;
    const amrex::Real q = getCharge();

    // Prepare interpolation of current components to cell center
    // The arrays below store the index type (staggering) of each MultiFab, with the third
    // component set to zero in the two-dimensional case.
    auto j_Nodal_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto jx_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto jy_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto jz_CC_type = amrex::GpuArray<int, 3>{0, 0, 0};
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        j_Nodal_type[i] = tmp_jx_fluid.ixType()[i];
        jx_CC_type[i] = jx.ixType()[i];
        jy_CC_type[i] = jy.ixType()[i];
        jz_CC_type[i] = jz.ixType()[i];
    }

    // We now need to create a mask to fix the double counting.
    WarpX &warpx = WarpX::GetInstance();
    const amrex::Geometry &geom = warpx.Geom(lev);
    const amrex::Periodicity &period = geom.periodicity();
    auto const &owner_mask_x = amrex::OwnerMask(jx, period);
    auto const &owner_mask_y = amrex::OwnerMask(jy, period);
    auto const &owner_mask_z = amrex::OwnerMask(jz, period);

// Calculate j at the nodes
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

        amrex::Array4<amrex::Real> tmp_jx_fluid_arr = tmp_jx_fluid.array(mfi);
        amrex::Array4<amrex::Real> tmp_jy_fluid_arr = tmp_jy_fluid.array(mfi);
        amrex::Array4<amrex::Real> tmp_jz_fluid_arr = tmp_jz_fluid.array(mfi);

        // Loop over cells (ParallelFor)
        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Calculate J from fluid and add it to jx/jy/jz
                auto gamma = std::sqrt(N_arr(i, j, k) * N_arr(i, j, k) + (NUx_arr(i, j, k) * NUx_arr(i, j, k) + NUy_arr(i, j, k) * NUy_arr(i, j, k) + NUz_arr(i, j, k) * NUz_arr(i, j, k)) * inv_clight_sq) / N_arr(i, j, k);
                tmp_jx_fluid_arr(i, j, k) = q * (NUx_arr(i, j, k) / gamma);
                tmp_jy_fluid_arr(i, j, k) = q * (NUy_arr(i, j, k) / gamma);
                tmp_jz_fluid_arr(i, j, k) = q * (NUz_arr(i, j, k) / gamma);
            }
        );
    }

// Interpolate j from the nodes to the simulation mesh (typically Yee mesh)
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box const &tile_box_cc_x = mfi.tilebox(jx.ixType().toIntVect());
        amrex::Box const &tile_box_cc_y = mfi.tilebox(jy.ixType().toIntVect());
        amrex::Box const &tile_box_cc_z = mfi.tilebox(jz.ixType().toIntVect());

        amrex::Array4<amrex::Real> jx_arr = jx.array(mfi);
        amrex::Array4<amrex::Real> jy_arr = jy.array(mfi);
        amrex::Array4<amrex::Real> jz_arr = jz.array(mfi);

        amrex::Array4<amrex::Real> tmp_jx_fluid_arr = tmp_jx_fluid.array(mfi);
        amrex::Array4<amrex::Real> tmp_jy_fluid_arr = tmp_jy_fluid.array(mfi);
        amrex::Array4<amrex::Real> tmp_jz_fluid_arr = tmp_jz_fluid.array(mfi);

        auto owner_mask_x_arr = owner_mask_x->array(mfi);
        auto owner_mask_y_arr = owner_mask_y->array(mfi);
        auto owner_mask_z_arr = owner_mask_z->array(mfi);

        // When using the `Interp` function, one needs to specify whether coarsening is desired.
        // Here, we do not perform any coarsening.
        amrex::GpuArray<int, 3U> coarsening_ratio = {AMREX_D_DECL(1, 1, 1)};

        // Loop over cells (ParallelFor)
        amrex::ParallelFor( tile_box_cc_x, tile_box_cc_y, tile_box_cc_z,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real jx_CC = ablastr::coarsen::sample::Interp(tmp_jx_fluid_arr,
                    j_Nodal_type, jx_CC_type, coarsening_ratio, i, j, k, 0);
                if ( owner_mask_x_arr(i,j,k) ) jx_arr(i, j, k) += jx_CC;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real jy_CC = ablastr::coarsen::sample::Interp(tmp_jy_fluid_arr,
                    j_Nodal_type, jy_CC_type, coarsening_ratio, i, j, k, 0);
                if ( owner_mask_y_arr(i,j,k) ) jy_arr(i, j, k) += jy_CC;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real jz_CC = ablastr::coarsen::sample::Interp(tmp_jz_fluid_arr,
                    j_Nodal_type, jz_CC_type, coarsening_ratio, i, j, k, 0);
                if ( owner_mask_z_arr(i,j,k) ) jz_arr(i, j, k) += jz_CC;
            }
        );
    }
}
