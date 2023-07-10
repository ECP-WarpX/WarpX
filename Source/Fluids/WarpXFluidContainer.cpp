/* Copyright 2023 Grant Johnson, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ablastr/coarsen/sample.H"
#include "Particles/Pusher/UpdateMomentumHigueraCary.H"
#include "Utils/WarpXProfilerWrapper.H"

#include "MusclHancockUtils.H"
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

    WarpX &warpx = WarpX::GetInstance();

    // set human-readable tag for each MultiFab
    auto const tag = [lev](std::string tagname)
    {
        tagname.append("[l=").append(std::to_string(lev)).append("]");
        return tagname;
    };

    warpx.AllocInitMultiFab(N[lev], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid density"), 0.0_rt);

    warpx.AllocInitMultiFab(NU[lev][0], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid momentum density [x]"), 0.0_rt);
    warpx.AllocInitMultiFab(NU[lev][1], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid momentum density [y]"), 0.0_rt);
    warpx.AllocInitMultiFab(NU[lev][2], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid momentum density [z]"), 0.0_rt);
}

void WarpXFluidContainer::InitData(int lev)
{
    WARPX_PROFILE("WarpXFluidContainer::InitData");

    // Extract objects that give the initial density and momentum
    InjectorDensity *inj_rho = plasma_injector->getInjectorDensity();
    InjectorMomentum *inj_mom = plasma_injector->getInjectorMomentumDevice();

    // Extract grid geometry properties
    WarpX &warpx = WarpX::GetInstance();
    const amrex::Geometry &geom = warpx.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const amrex::Real clight = PhysConst::c;

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
                NUx_arr(i, j, k) = n * u.x * clight;
                NUy_arr(i, j, k) = n * u.y * clight;
                NUz_arr(i, j, k) = n * u.z * clight;
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

    // Step the Lorentz Term
    GatherAndPush(lev, Ex, Ey, Ez, Bx, By, Bz);

    // Step the Advective term
    AdvectivePush_Muscl(lev);

    // Deposit J to the simulation mesh
    if (!skip_deposition)
    {
        DepositCurrent(lev, jx, jy, jz);
    }
}


// Muscl Advection Update
void WarpXFluidContainer::AdvectivePush_Muscl (int lev)
{
    WARPX_PROFILE("WarpXFluidContainer::AdvectivePush_Muscl");

    // Grab the grid spacing
    WarpX &warpx = WarpX::GetInstance();
    const Real dt = warpx.getdt(lev);
    const amrex::Geometry &geom = warpx.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const amrex::Real clight = PhysConst::c;
    const amrex::Periodicity &period = geom.periodicity();
    auto cx = (dt/dx[0]);
    auto cy = (dt/dx[1]);
    auto cz = (dt/dx[2]);
    auto cx_half = 0.5*(dt/dx[0]);
    auto cy_half = 0.5*(dt/dx[1]);
    auto cz_half = 0.5*(dt/dx[2]);

    // Advection push
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

                // - Grab local Uz Uy Ux gamma
                // Isolate U from NU
                auto Ux = (NUx_arr(i, j, k) / N_arr(i,j,k));
                auto Uy = (NUy_arr(i, j, k) / N_arr(i,j,k));
                auto Uz = (NUz_arr(i, j, k) / N_arr(i,j,k));
                auto gamma = sqrt(1.0 + (Ux*Ux + Uy*Uy + Uz*Uz)/(clight*clight) );
                auto a = amrex::Math::powi<2>( clight * gamma ) * gamma;

                // Grab boundaries x
                auto Ux_Rx = (NUx_arr(i+1, j, k) / N_arr(i+1,j,k));
                auto Uy_Rx = (NUy_arr(i+1, j, k) / N_arr(i+1,j,k));
                auto Uz_Rx = (NUz_arr(i+1, j, k) / N_arr(i+1,j,k));
                auto gamma_Rx = sqrt(1.0 + (Ux_Rx*Ux_Rx + Uy_Rx*Uy_Rx + Uz_Rx*Uz_Rx)/(clight*clight) );
                auto a_Rx = amrex::Math::powi<2>( clight * gamma_Rx ) * gamma_Rx;
                auto Ux_Lx = (NUx_arr(i-1, j, k) / N_arr(i-1,j,k));
                auto Uy_Lx = (NUy_arr(i-1, j, k) / N_arr(i-1,j,k));
                auto Uz_Lx = (NUz_arr(i-1, j, k) / N_arr(i-1,j,k));
                auto gamma_Lx = sqrt(1.0 + (Ux_Lx*Ux_Lx + Uy_Lx*Uy_Lx + Uz_Lx*Uz_Lx)/(clight*clight) );
                auto a_Lx = amrex::Math::powi<2>( clight * gamma_Lx ) * gamma_Lx;

                // Grab boundaries y
                auto Ux_Ry = (NUx_arr(i,j+1,k) / N_arr(i,j+1,k));
                auto Uy_Ry = (NUy_arr(i,j+1,k) / N_arr(i,j+1,k));
                auto Uz_Ry = (NUz_arr(i,j+1,k) / N_arr(i,j+1,k));
                auto gamma_Ry = sqrt(1.0 + (Ux_Ry*Ux_Ry + Uy_Ry*Uy_Ry + Uz_Ry*Uz_Ry)/(clight*clight) );
                auto a_Ry = amrex::Math::powi<2>( clight * gamma_Ry ) * gamma_Ry;
                auto Ux_Ly = (NUx_arr(i,j-1,k) / N_arr(i,j-1,k));
                auto Uy_Ly = (NUy_arr(i,j-1,k) / N_arr(i,j-1,k));
                auto Uz_Ly = (NUz_arr(i,j-1,k) / N_arr(i,j-1,k));
                auto gamma_Ly = sqrt(1.0 + (Ux_Ly*Ux_Ly + Uy_Ly*Uy_Ly + Uz_Ly*Uz_Ly)/(clight*clight) );
                auto a_Ly = amrex::Math::powi<2>( clight * gamma_Ly ) * gamma_Ly;

                // Grab boundaries z
                auto Ux_Rz = (NUx_arr(i,j,k+1) / N_arr(i,j,k+1));
                auto Uy_Rz = (NUy_arr(i,j,k+1) / N_arr(i,j,k+1));
                auto Uz_Rz = (NUz_arr(i,j,k+1) / N_arr(i,j,k+1));
                auto gamma_Rz = sqrt(1.0 + (Ux_Rz*Ux_Rz + Uy_Rz*Uy_Rz + Uz_Rz*Uz_Rz)/(clight*clight) );
                auto a_Rz = amrex::Math::powi<2>( clight * gamma_Rz ) * gamma_Rz;
                auto Ux_Lz = (NUx_arr(i,j,k-1) / N_arr(i,j,k-1));
                auto Uy_Lz = (NUy_arr(i,j,k-1) / N_arr(i,j,k-1));
                auto Uz_Lz = (NUz_arr(i,j,k-1) / N_arr(i,j,k-1));
                auto gamma_Lz = sqrt(1.0 + (Ux_Lz*Ux_Lz + Uy_Lz*Uy_Lz + Uz_Lz*Uz_Lz)/(clight*clight) );
                auto a_Lz = amrex::Math::powi<2>( clight * gamma_Lz ) * gamma_Lz;

                // Memory for midpoint values
                amrex::Real Q_minus_x[4];
                amrex::Real Q_minus_y[4];
                amrex::Real Q_minus_z[4];
                amrex::Real Q_plus_x[4];
                amrex::Real Q_plus_y[4];
                amrex::Real Q_plus_z[4];

                amrex::Real Q_minus_xL[4];
                amrex::Real Q_minus_yL[4];
                amrex::Real Q_minus_zL[4];
                amrex::Real Q_plus_xR[4];
                amrex::Real Q_plus_yR[4];
                amrex::Real Q_plus_zR[4];

                amrex::Real tmp_unused[4];

                // Select the specific implmentation depending on dimensionality
                #if defined(WARPX_DIM_3D)

                // Grab Q's required
                Q_step(i, j, k, cx_half, cy_half, cz_half, Ux, Uy, Uz, a,
                clight, gamma, N_arr,  NUx_arr, NUy_arr, NUz_arr,
                Q_minus_x, Q_minus_y, Q_minus_z, Q_plus_x, Q_plus_y, Q_plus_z);

                Q_step(i-1, j, k, cx_half, cy_half, cz_half, Ux_Lx, Uy_Lx, Uz_Lx, a_Lx,
                clight, gamma_Lx, N_arr,  NUx_arr, NUy_arr, NUz_arr,
                Q_minus_xL, tmp_unused, tmp_unused, tmp_unused, tmp_unused, tmp_unused);
                Q_step(i+1, j, k, cx_half, cy_half, cz_half, Ux_Rx, Uy_Rx, Uz_Rx, a_Rx,
                clight, gamma_Rx, N_arr,  NUx_arr, NUy_arr, NUz_arr,
                tmp_unused, tmp_unused, tmp_unused, Q_plus_xR, tmp_unused, tmp_unused);

                Q_step(i, j-1, k, cx_half, cy_half, cz_half, Ux_Ly, Uy_Ly, Uz_Ly, a_Ly,
                clight, gamma_Ly, N_arr,  NUx_arr, NUy_arr, NUz_arr,
                tmp_unused, Q_minus_yL, tmp_unused, tmp_unused, tmp_unused, tmp_unused);
                Q_step(i, j+1, k, cx_half, cy_half, cz_half, Ux_Ry, Uy_Ry, Uz_Ry, a_Ry,
                clight, gamma_Ry, N_arr,  NUx_arr, NUy_arr, NUz_arr,
                tmp_unused, tmp_unused, tmp_unused, tmp_unused, Q_plus_yR, tmp_unused);

                Q_step(i, j, k-1, cx_half, cy_half, cz_half, Ux_Lz, Uy_Lz, Uz_Lz, a_Lz,
                clight, gamma_Lz, N_arr,  NUx_arr, NUy_arr, NUz_arr,
                tmp_unused, tmp_unused, Q_minus_zL, tmp_unused, tmp_unused, tmp_unused);
                Q_step(i, j, k+1, cx_half, cy_half, cz_half, Ux_Rz, Uy_Rz, Uz_Rz, a_Rz,
                clight, gamma_Rz, N_arr,  NUx_arr, NUy_arr, NUz_arr,
                tmp_unused, tmp_unused, tmp_unused, tmp_unused, tmp_unused, Q_plus_zR);

                // Compute Vx Vy Vz (R is right (+1), L is left (-1))
                auto Vx = Ux/gamma;
                auto Vy = Uy/gamma;
                auto Vz = Uz/gamma;
                auto Vx_R = Ux_Rx/gamma_Rx;
                auto Vy_R = Uy_Ry/gamma_Ry;
                auto Vz_R = Uz_Rz/gamma_Rz;
                auto Vx_L = Ux_Lx/gamma_Lx;
                auto Vy_L = Uy_Ly/gamma_Ly;
                auto Vz_L = Uz_Lz/gamma_Lz;

                // compute the fluxes:
                auto F0_minusx = flux(Q_minus_xL[0],Q_plus_x[0],  Vx_L, Vx);
                auto F0_plusx =  flux(Q_minus_x[0], Q_plus_xR[0],Vx,   Vx_R);
                auto F1_minusx = flux(Q_minus_xL[1],Q_plus_x[1],  Vx_L, Vx);
                auto F1_plusx =  flux(Q_minus_x[1], Q_plus_xR[1],Vx,   Vx_R);
                auto F2_minusx = flux(Q_minus_xL[2],Q_plus_x[2],  Vx_L, Vx);
                auto F2_plusx =  flux(Q_minus_x[2], Q_plus_xR[2],Vx,   Vx_R);
                auto F3_minusx = flux(Q_minus_xL[3],Q_plus_x[3],  Vx_L, Vx);
                auto F3_plusx =  flux(Q_minus_x[3], Q_plus_xR[3],Vx,   Vx_R);

                auto F0_minusy = flux(Q_minus_yL[0],Q_plus_y[0],  Vy_L,Vy);
                auto F0_plusy =  flux(Q_minus_y[0], Q_plus_yR[0],Vy,  Vy_R);
                auto F1_minusy = flux(Q_minus_yL[1],Q_plus_y[1],  Vy_L,Vy);
                auto F1_plusy =  flux(Q_minus_y[1], Q_plus_yR[1],Vy,  Vy_R);
                auto F2_minusy = flux(Q_minus_yL[2],Q_plus_y[2],  Vy_L,Vy);
                auto F2_plusy =  flux(Q_minus_y[2], Q_plus_yR[2],Vy,  Vy_R);
                auto F3_minusy = flux(Q_minus_yL[3],Q_plus_y[3],  Vy_L,Vy);
                auto F3_plusy =  flux(Q_minus_y[3], Q_plus_yR[3],Vy,  Vy_R);

                auto F0_minusz = flux(Q_minus_zL[0],Q_plus_z[0],  Vz_L,Vz);
                auto F0_plusz =  flux(Q_minus_z[0], Q_plus_zR[0],Vz,  Vz_R);
                auto F1_minusz = flux(Q_minus_zL[1],Q_plus_z[1],  Vz_L,Vz);
                auto F1_plusz =  flux(Q_minus_z[1], Q_plus_zR[1],Vz,  Vz_R);
                auto F2_minusz = flux(Q_minus_zL[2],Q_plus_z[2],  Vz_L,Vz);
                auto F2_plusz =  flux(Q_minus_z[2], Q_plus_zR[2],Vz,  Vz_R);
                auto F3_minusz = flux(Q_minus_zL[3],Q_plus_z[3],  Vz_L,Vz);
                auto F3_plusz =  flux(Q_minus_z[3], Q_plus_zR[3],Vz,  Vz_R);

                // Update Q from tn -> tn + dt
                N_arr(i,j,k) = N_arr(i,j,k)     - cx*(F0_plusx - F0_minusx)
                                                - cy*(F0_plusy - F0_minusy)
                                                - cz*(F0_plusz - F0_minusz);
                NUx_arr(i,j,k) = NUx_arr(i,j,k) - cx*(F1_plusx - F1_minusx)
                                                - cy*(F1_plusy - F1_minusy)
                                                - cz*(F1_plusz - F1_minusz);
                NUy_arr(i,j,k) = NUy_arr(i,j,k) - cx*(F2_plusx - F2_minusx)
                                                - cy*(F2_plusy - F2_minusy)
                                                - cz*(F2_plusz - F2_minusz);
                NUz_arr(i,j,k) = NUz_arr(i,j,k) - cx*(F3_plusx - F3_minusx)
                                                - cy*(F3_plusy - F3_minusy)
                                                - cz*(F3_plusz - F3_minusz);

                #elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)

                #else

                #endif
            }
        );
    }

    // Fill guard cells
    FillBoundary(*N[lev], N[lev]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][0], NU[lev][0]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][1], NU[lev][1]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][2], NU[lev][2]->nGrowVect(), WarpX::do_single_precision_comms, period);
}

// Momentum source from fields
void WarpXFluidContainer::GatherAndPush (
    int lev,
    const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
    const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz)
{
    WARPX_PROFILE("WarpXFluidContainer::GatherAndPush");

    WarpX &warpx = WarpX::GetInstance();
    const amrex::Real q = getCharge();
    const amrex::Real m = getMass();
    const Real dt = warpx.getdt(lev);

   // Prepare interpolation of current components to cell center
    auto Nodal_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Ex_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Ey_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Ez_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Bx_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto By_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto Bz_type = amrex::GpuArray<int, 3>{0, 0, 0};
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        Nodal_type[i] = N[lev]->ixType()[i];
        Ex_type[i] = Ex.ixType()[i];
        Ey_type[i] = Ey.ixType()[i];
        Ez_type[i] = Ez.ixType()[i];
        Bx_type[i] = Bx.ixType()[i];
        By_type[i] = By.ixType()[i];
        Bz_type[i] = Bz.ixType()[i];
    }


    // H&C push the momentum
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

        amrex::Array4<const amrex::Real> const& Ex_arr = Ex.array(mfi);
        amrex::Array4<const amrex::Real> const& Ey_arr = Ey.array(mfi);
        amrex::Array4<const amrex::Real> const& Ez_arr = Ez.array(mfi);
        amrex::Array4<const amrex::Real> const& Bx_arr = Bx.array(mfi);
        amrex::Array4<const amrex::Real> const& By_arr = By.array(mfi);
        amrex::Array4<const amrex::Real> const& Bz_arr = Bz.array(mfi);

        // Here, we do not perform any coarsening.
        amrex::GpuArray<int, 3U> coarsening_ratio = {AMREX_D_DECL(1, 1, 1)};

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // Interpolate fields from tmp to Nodal points
                amrex::Real Ex_Nodal = ablastr::coarsen::sample::Interp(Ex_arr,
                    Ex_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Ey_Nodal = ablastr::coarsen::sample::Interp(Ey_arr,
                    Ey_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Ez_Nodal = ablastr::coarsen::sample::Interp(Ez_arr,
                    Ez_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Bx_Nodal = ablastr::coarsen::sample::Interp(Bx_arr,
                    Bx_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real By_Nodal = ablastr::coarsen::sample::Interp(By_arr,
                    By_type, Nodal_type, coarsening_ratio, i, j, k, 0);
                amrex::Real Bz_Nodal = ablastr::coarsen::sample::Interp(Bz_arr,
                    Bz_type, Nodal_type, coarsening_ratio, i, j, k, 0);

                // Isolate U from NU
                auto tmp_Ux = (NUx_arr(i, j, k) / N_arr(i,j,k));
                auto tmp_Uy = (NUy_arr(i, j, k) / N_arr(i,j,k));
                auto tmp_Uz = (NUz_arr(i, j, k) / N_arr(i,j,k));

                // Push the fluid momentum
                UpdateMomentumHigueraCary(tmp_Ux, tmp_Uy, tmp_Uz,
                    Ex_Nodal, Ey_Nodal, Ez_Nodal,
                    Bx_Nodal, By_Nodal, Bz_Nodal, q, m, dt );

                // Calculate NU
                NUx_arr(i,j,k) = N_arr(i,j,k)*tmp_Ux;
                NUy_arr(i,j,k) = N_arr(i,j,k)*tmp_Uy;
                NUz_arr(i,j,k) = N_arr(i,j,k)*tmp_Uz;

            }
        );
    }

    // Fill guard cells
    const amrex::Geometry &geom = warpx.Geom(lev);
    const amrex::Periodicity &period = geom.periodicity();
    FillBoundary(*NU[lev][0], NU[lev][0]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][1], NU[lev][1]->nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(*NU[lev][2], NU[lev][2]->nGrowVect(), WarpX::do_single_precision_comms, period);
}

void WarpXFluidContainer::DepositCharge(int lev, amrex::MultiFab &rho)
{
    WARPX_PROFILE("WarpXFluidContainer::DepositCharge");

    WarpX &warpx = WarpX::GetInstance();
    const amrex::Geometry &geom = warpx.Geom(lev);
    const amrex::Periodicity &period = geom.periodicity();
    const amrex::Real q = getCharge();
    auto const &owner_mask_rho = amrex::OwnerMask(rho, period);

    // Assertion, make sure rho is at the same location as N
    AMREX_ALWAYS_ASSERT(rho.ixType().nodeCentered());

    // Loop over and deposit charge density
    #ifdef AMREX_USE_OMP
    #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    #endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        amrex::Box const &tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());
        amrex::Array4<Real> const &N_arr = N[lev]->array(mfi);
        amrex::Array4<amrex::Real> rho_arr = rho.array(mfi);
        auto owner_mask_rho_arr = owner_mask_rho->array(mfi);

        // Deposit Rho
        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if ( owner_mask_rho_arr(i,j,k) ) rho_arr(i,j,k) += q*N_arr(i,j,k);
            }
        );
    }
}


void WarpXFluidContainer::DepositCurrent(
    int lev,
    amrex::MultiFab &jx, amrex::MultiFab &jy, amrex::MultiFab &jz)
{
    WARPX_PROFILE("WarpXFluidContainer::DepositCurrent");

    // Temporary nodal currents
    amrex::MultiFab tmp_jx_fluid(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_jy_fluid(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);
    amrex::MultiFab tmp_jz_fluid(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 0);

    const amrex::Real inv_clight_sq = 1.0_prt / PhysConst::c / PhysConst::c;
    const amrex::Real q = getCharge();

    // Prepare interpolation of current components to cell center
    auto j_nodal_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto jx_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto jy_type = amrex::GpuArray<int, 3>{0, 0, 0};
    auto jz_type = amrex::GpuArray<int, 3>{0, 0, 0};
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        j_nodal_type[i] = tmp_jx_fluid.ixType()[i];
        jx_type[i] = jx.ixType()[i];
        jy_type[i] = jy.ixType()[i];
        jz_type[i] = jz.ixType()[i];
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

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Calculate J from fluid quantities
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
        amrex::Box const &tile_box_x = mfi.tilebox(jx.ixType().toIntVect());
        amrex::Box const &tile_box_y = mfi.tilebox(jy.ixType().toIntVect());
        amrex::Box const &tile_box_z = mfi.tilebox(jz.ixType().toIntVect());

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

        // Interpolate fluid current and deposit it
        // ( mask double counting )
        amrex::ParallelFor( tile_box_x, tile_box_y, tile_box_z,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real jx_tmp = ablastr::coarsen::sample::Interp(tmp_jx_fluid_arr,
                    j_nodal_type, jx_type, coarsening_ratio, i, j, k, 0);
                if ( owner_mask_x_arr(i,j,k) ) jx_arr(i, j, k) += jx_tmp;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real jy_tmp = ablastr::coarsen::sample::Interp(tmp_jy_fluid_arr,
                    j_nodal_type, jy_type, coarsening_ratio, i, j, k, 0);
                if ( owner_mask_y_arr(i,j,k) ) jy_arr(i, j, k) += jy_tmp;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real jz_tmp = ablastr::coarsen::sample::Interp(tmp_jz_fluid_arr,
                    j_nodal_type, jz_type, coarsening_ratio, i, j, k, 0);
                if ( owner_mask_z_arr(i,j,k) ) jz_arr(i, j, k) += jz_tmp;
            }
        );
    }
}
