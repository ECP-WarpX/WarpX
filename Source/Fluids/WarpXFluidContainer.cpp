/* Copyright 2023 Grant Johnson, Remi Lehe
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


    /** Todo */
    // - Grab local Uz Uy Ux gamma
    // - Compute the Flux Jacobian (A = dF/dQ) in each direction
    // - Compute the cell slopes
        // - ave function required ( minmod )
    // - Predict Q at dt/2
    // - Predict Q at cell edges using the cell slopes
    // - Update Q at t + dt
    // - Update N, NU

    // Temporary velocities
    amrex::MultiFab tmp_Vx(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Vy(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Vz(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);

    // Temporary Half-step values
    amrex::MultiFab tmp_Q_minus1_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus2_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus3_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus4_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus1_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus2_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus3_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus4_x(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus1_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus2_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus3_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus4_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus1_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus2_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus3_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus4_y(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus1_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus2_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus3_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_minus4_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus1_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus2_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus3_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);
    amrex::MultiFab tmp_Q_plus4_z(N[lev]->boxArray(), N[lev]->DistributionMap(), 1, 1);

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

        amrex::Array4<amrex::Real> Vx = tmp_Vx.array(mfi);
        amrex::Array4<amrex::Real> Vy = tmp_Vy.array(mfi);
        amrex::Array4<amrex::Real> Vz = tmp_Vz.array(mfi);

        amrex::Array4<amrex::Real> Q_minus1_x = tmp_Q_minus1_x.array(mfi);
        amrex::Array4<amrex::Real> Q_minus2_x = tmp_Q_minus2_x.array(mfi);
        amrex::Array4<amrex::Real> Q_minus3_x = tmp_Q_minus3_x.array(mfi);
        amrex::Array4<amrex::Real> Q_minus4_x = tmp_Q_minus4_x.array(mfi);
        amrex::Array4<amrex::Real> Q_plus1_x = tmp_Q_plus1_x.array(mfi);
        amrex::Array4<amrex::Real> Q_plus2_x = tmp_Q_plus2_x.array(mfi);
        amrex::Array4<amrex::Real> Q_plus3_x = tmp_Q_plus3_x.array(mfi);
        amrex::Array4<amrex::Real> Q_plus4_x = tmp_Q_plus4_x.array(mfi);

        amrex::Array4<amrex::Real> Q_minus1_y = tmp_Q_minus1_y.array(mfi);
        amrex::Array4<amrex::Real> Q_minus2_y = tmp_Q_minus2_y.array(mfi);
        amrex::Array4<amrex::Real> Q_minus3_y = tmp_Q_minus3_y.array(mfi);
        amrex::Array4<amrex::Real> Q_minus4_y = tmp_Q_minus4_y.array(mfi);
        amrex::Array4<amrex::Real> Q_plus1_y = tmp_Q_plus1_y.array(mfi);
        amrex::Array4<amrex::Real> Q_plus2_y = tmp_Q_plus2_y.array(mfi);
        amrex::Array4<amrex::Real> Q_plus3_y = tmp_Q_plus3_y.array(mfi);
        amrex::Array4<amrex::Real> Q_plus4_y = tmp_Q_plus4_y.array(mfi);

        amrex::Array4<amrex::Real> Q_minus1_z = tmp_Q_minus1_z.array(mfi);
        amrex::Array4<amrex::Real> Q_minus2_z = tmp_Q_minus2_z.array(mfi);
        amrex::Array4<amrex::Real> Q_minus3_z = tmp_Q_minus3_z.array(mfi);
        amrex::Array4<amrex::Real> Q_minus4_z = tmp_Q_minus4_z.array(mfi);
        amrex::Array4<amrex::Real> Q_plus1_z = tmp_Q_plus1_z.array(mfi);
        amrex::Array4<amrex::Real> Q_plus2_z = tmp_Q_plus2_z.array(mfi);
        amrex::Array4<amrex::Real> Q_plus3_z = tmp_Q_plus3_z.array(mfi);
        amrex::Array4<amrex::Real> Q_plus4_z = tmp_Q_plus4_z.array(mfi);


        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // - Grab local Uz Uy Ux gamma
                // Isolate U from NU
                auto Ux = (NUx_arr(i, j, k) / N_arr(i,j,k));
                auto Uy = (NUy_arr(i, j, k) / N_arr(i,j,k));
                auto Uz = (NUz_arr(i, j, k) / N_arr(i,j,k));
                auto gamma = sqrt(1.0 + (Ux*Ux + Uy*Uy + Uz*Uz)/(clight*clight) );
                auto a = (clight*clight*clight*clight)+(2.0*(Uz*Uz)+2.0*(Uy*Uy)+
                2.0*(Ux*Ux))*(clight*clight)+(Uz*Uz*Uz*Uz)+(2.0*(Uy*Uy)+2.0*(Ux*Ux))
                *(Uz*Uz)+(Uy*Uy*Uy*Uy)+2.0*(Ux*Ux)*(Uy*Uy)+(Ux*Ux*Ux*Ux);

                // Compute Vx Vy Vz
                Vx(i,j,k) = Ux/gamma;
                Vy(i,j,k) = Uy/gamma;
                Vz(i,j,k) = Uz/gamma;

                // Select the specific implmentation depending on dimensionality
                #if defined(WARPX_DIM_3D)

                // Compute the Flux-Jacobian Elements in x
                auto A11x = ((Ux*(Uz*Uz)+Ux*(Uy*Uy)+(Ux*Ux*Ux))*(clight*clight)*gamma)/a;
                auto A12x = (((clight*clight*clight*clight)+((Uz*Uz)+(Uy*Uy))*(clight*clight))*gamma)/a;
                auto A13x = -(Ux*Uy)/((clight*clight)*(gamma*gamma*gamma));
                auto A14x = -(Ux*Uz)/((clight*clight)*(gamma*gamma*gamma));

                auto A21x = -(Ux*Ux)/(gamma*gamma*gamma);
                auto A22x = ((2.0*Ux*(clight*clight*clight*clight)+(2.0*Ux*(Uz*Uz)+2.0*Ux*(Uy*Uy)+(Ux*Ux*Ux))*(clight*clight))*gamma)/a;
                auto A23x = -((Ux*Ux)*Uy)/((clight*clight)*(gamma*gamma*gamma));
                auto A24x = -((Ux*Ux)*Uz)/((clight*clight)*(gamma*gamma*gamma));

                auto A31x = -(Ux*Uy)/(gamma*gamma*gamma);
                auto A32x = ((Uy*(clight*clight*clight*clight)+(Uy*(Uz*Uz)+(Uy*Uy*Uy))*(clight*clight))*gamma)/a;
                auto A33x = ((Ux*(clight*clight*clight*clight)+(Ux*(Uz*Uz)+(Ux*Ux*Ux))*(clight*clight))*gamma)/a;
                auto A34x = -(Ux*Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));

                auto A41x = -(Ux*Uz)/(gamma*gamma*gamma);
                auto A42x = ((Uz*(clight*clight*clight*clight)+((Uz*Uz*Uz)+(Uy*Uy)*Uz)*(clight*clight))*gamma)/a;
                auto A43x = -(Ux*Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));
                auto A44x = ((Ux*(clight*clight*clight*clight)+(Ux*(Uy*Uy)+(Ux*Ux*Ux))*(clight*clight))*gamma)/a;


                // Compute the Flux-Jacobian Elements in y
                auto A11y = ((Uy*(Uz*Uz)+(Uy*Uy*Uy)+(Ux*Ux)*Uy)*(clight*clight)*gamma)/a;
                auto A12y = -(Ux*Uy)/((clight*clight)*(gamma*gamma*gamma));
                auto A13y = (((clight*clight*clight*clight)+((Uz*Uz)+(Ux*Ux))*(clight*clight))*gamma)/a;
                auto A14y = -(Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));

                auto A21y = -(Ux*Uy)/(gamma*gamma*gamma);
                auto A22y = ((Uy*(clight*clight*clight*clight)+(Uy*(Uz*Uz)+(Uy*Uy*Uy))*(clight*clight))*gamma)/a;
                auto A23y = ((Ux*(clight*clight*clight*clight)+(Ux*(Uz*Uz)+(Ux*Ux*Ux))*(clight*clight))*gamma)/a;
                auto A24y = -(Ux*Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));

                auto A31y = -(Uy*Uy)/(gamma*gamma*gamma);
                auto A32y = -(Ux*(Uy*Uy))/((clight*clight)*(gamma*gamma*gamma));
                auto A33y = ((2.0*Uy*(clight*clight*clight*clight)+(2.0*Uy*(Uz*Uz)+(Uy*Uy*Uy)+2.0*(Ux*Ux)*Uy)*(clight*clight))*gamma)/a;
                auto A34y = -((Uy*Uy)*Uz)/((clight*clight)*(gamma*gamma*gamma));

                auto A41y = -(Uy*Uz)/(gamma*gamma*gamma);
                auto A42y = -(Ux*Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));
                auto A43y = ((Uz*(clight*clight*clight*clight)+((Uz*Uz*Uz)+(Ux*Ux)*Uz)*(clight*clight))*gamma)/a;
                auto A44y = ((Uy*(clight*clight*clight*clight)+((Uy*Uy*Uy)+(Ux*Ux)*Uy)*(clight*clight))*gamma)/a;


                // Compute the Flux-Jacobian Elements in z
                auto A11z = (((Uz*Uz*Uz)+((Uy*Uy)+(Ux*Ux))*Uz)*(clight*clight)*gamma)/a;
                auto A12z = -(Ux*Uz)/((clight*clight)*(gamma*gamma*gamma));
                auto A13z = -(Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));
                auto A14z = (((clight*clight*clight*clight)+((Uy*Uy)+(Ux*Ux))*(clight*clight))*gamma)/a;

                auto A21z = -(Ux*Uz)/(gamma*gamma*gamma);
                auto A22z = ((Uz*(clight*clight*clight*clight)+((Uz*Uz*Uz)+(Uy*Uy)*Uz)*(clight*clight))*gamma)/a;
                auto A23z = -(Ux*Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));
                auto A24z = ((Ux*(clight*clight*clight*clight)+(Ux*(Uy*Uy)+(Ux*Ux*Ux))*(clight*clight))*gamma)/a;

                auto A31z = -(Uy*Uz)/(gamma*gamma*gamma);
                auto A32z = -(Ux*Uy*Uz)/((clight*clight)*(gamma*gamma*gamma));
                auto A33z = ((Uz*(clight*clight*clight*clight)+((Uz*Uz*Uz)+(Ux*Ux)*Uz)*(clight*clight))*gamma)/a;
                auto A34z = ((Uy*(clight*clight*clight*clight)+((Uy*Uy*Uy)+(Ux*Ux)*Uy)*(clight*clight))*gamma)/a;

                auto A41z = -(Uz*Uz)/(gamma*gamma*gamma);
                auto A42z = -(Ux*(Uz*Uz))/((clight*clight)*(gamma*gamma*gamma));
                auto A43z = -(Uy*(Uz*Uz))/((clight*clight)*(gamma*gamma*gamma));
                auto A44z = ((2.0*Uz*(clight*clight*clight*clight)+((Uz*Uz*Uz)+(2.0*(Uy*Uy)+2.0*(Ux*Ux))*Uz)*(clight*clight))*gamma)/a;

                // Compute the cell slopes x
                auto dQ1x = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                auto dQ2x = ave( NUx_arr(i,j,k) - NUx_arr(i-1,j,k) , NUx_arr(i+1,j,k) - NUx_arr(i,j,k) );
                auto dQ3x = ave( NUy_arr(i,j,k) - NUy_arr(i-1,j,k) , NUy_arr(i+1,j,k) - NUy_arr(i,j,k) );
                auto dQ4x = ave( NUz_arr(i,j,k) - NUz_arr(i-1,j,k) , NUz_arr(i+1,j,k) - NUz_arr(i,j,k) );

                // Compute the cell slopes y
                auto dQ1y = ave( N_arr(i,j,k) - N_arr(i,j-1,k) , N_arr(i,j+1,k) - N_arr(i,j,k) );
                auto dQ2y = ave( NUx_arr(i,j,k) - NUx_arr(i,j-1,k) , NUx_arr(i,j+1,k) - NUx_arr(i,j,k) );
                auto dQ3y = ave( NUy_arr(i,j,k) - NUy_arr(i,j-1,k) , NUy_arr(i,j+1,k) - NUy_arr(i,j,k) );
                auto dQ4y = ave( NUz_arr(i,j,k) - NUz_arr(i,j-1,k) , NUz_arr(i,j+1,k) - NUz_arr(i,j,k) );

                // Compute the cell slopes z
                auto dQ1z = ave( N_arr(i,j,k) - N_arr(i,j,k-1) , N_arr(i,j,k+1) - N_arr(i,j,k) );
                auto dQ2z = ave( NUx_arr(i,j,k) - NUx_arr(i,j,k-1) , NUx_arr(i,j,k+1) - NUx_arr(i,j,k) );
                auto dQ3z = ave( NUy_arr(i,j,k) - NUy_arr(i,j,k-1) , NUy_arr(i,j,k+1) - NUy_arr(i,j,k) );
                auto dQ4z = ave( NUz_arr(i,j,k) - NUz_arr(i,j,k-1) , NUz_arr(i,j,k+1) - NUz_arr(i,j,k) );

                // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
                auto AdQ1x = A11x*dQ1x + A12x*dQ2x + A13x*dQ3x + A14x*dQ4x;
                auto AdQ2x = A21x*dQ1x + A22x*dQ2x + A23x*dQ3x + A24x*dQ4x;
                auto AdQ3x = A31x*dQ1x + A32x*dQ2x + A33x*dQ3x + A34x*dQ4x;
                auto AdQ4x = A41x*dQ1x + A42x*dQ2x + A43x*dQ3x + A44x*dQ4x;
                auto AdQ1y = A11y*dQ1y + A12y*dQ2y + A13y*dQ3y + A14y*dQ4y;
                auto AdQ2y = A21y*dQ1y + A22y*dQ2y + A23y*dQ3y + A24y*dQ4y;
                auto AdQ3y = A31y*dQ1y + A32y*dQ2y + A33y*dQ3y + A34y*dQ4y;
                auto AdQ4y = A41y*dQ1y + A42y*dQ2y + A43y*dQ3y + A44y*dQ4y;
                auto AdQ1z = A11z*dQ1z + A12z*dQ2z + A13z*dQ3z + A14z*dQ4z;
                auto AdQ2z = A21z*dQ1z + A22z*dQ2z + A23z*dQ3z + A24z*dQ4z;
                auto AdQ3z = A31z*dQ1z + A32z*dQ2z + A33z*dQ3z + A34z*dQ4z;
                auto AdQ4z = A41z*dQ1z + A42z*dQ2z + A43z*dQ3z + A44z*dQ4z;
                auto Q_tilde1 = N_arr(i,j,k)          - cx_half*AdQ1x - cy_half*AdQ1y - cz_half*AdQ1z;
                auto Q_tilde2 = NUx_arr(i,j,k) - cx_half*AdQ2x - cy_half*AdQ2y - cz_half*AdQ2z;
                auto Q_tilde3 = NUy_arr(i,j,k) - cx_half*AdQ3x - cy_half*AdQ3y - cz_half*AdQ3z;
                auto Q_tilde4 = NUz_arr(i,j,k) - cx_half*AdQ4x - cy_half*AdQ4y - cz_half*AdQ4z;

                // Predict Q at the cell edges (x)
                Q_minus1_x(i,j,k) = Q_tilde1 + dQ1x/2.0;
                Q_minus2_x(i,j,k) = Q_tilde2 + dQ2x/2.0;
                Q_minus3_x(i,j,k) = Q_tilde3 + dQ3x/2.0;
                Q_minus4_x(i,j,k) = Q_tilde4 + dQ4x/2.0;
                Q_plus1_x(i,j,k) = Q_tilde1 - dQ1x/2.0;
                Q_plus2_x(i,j,k) = Q_tilde2 - dQ2x/2.0;
                Q_plus3_x(i,j,k) = Q_tilde3 - dQ3x/2.0;
                Q_plus4_x(i,j,k) = Q_tilde4 - dQ4x/2.0;

                // Predict Q at the cell edges (y)
                Q_minus1_y(i,j,k) = Q_tilde1 + dQ1y/2.0;
                Q_minus2_y(i,j,k) = Q_tilde2 + dQ2y/2.0;
                Q_minus3_y(i,j,k) = Q_tilde3 + dQ3y/2.0;
                Q_minus4_y(i,j,k) = Q_tilde4 + dQ4y/2.0;
                Q_plus1_y(i,j,k) = Q_tilde1 - dQ1y/2.0;
                Q_plus2_y(i,j,k) = Q_tilde2 - dQ2y/2.0;
                Q_plus3_y(i,j,k) = Q_tilde3 - dQ3y/2.0;
                Q_plus4_y(i,j,k) = Q_tilde4 - dQ4y/2.0;

                // Predict Q at the cell edges (z)
                Q_minus1_z(i,j,k) = Q_tilde1 + dQ1z/2.0;
                Q_minus2_z(i,j,k) = Q_tilde2 + dQ2z/2.0;
                Q_minus3_z(i,j,k) = Q_tilde3 + dQ3z/2.0;
                Q_minus4_z(i,j,k) = Q_tilde4 + dQ4z/2.0;
                Q_plus1_z(i,j,k) = Q_tilde1 - dQ1z/2.0;
                Q_plus2_z(i,j,k) = Q_tilde2 - dQ2z/2.0;
                Q_plus3_z(i,j,k) = Q_tilde3 - dQ3z/2.0;
                Q_plus4_z(i,j,k) = Q_tilde4 - dQ4z/2.0;


                #elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)

                #else

                #endif
            }
        );
    }

    FillBoundary(tmp_Vx, tmp_Vx.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Vy, tmp_Vy.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Vz, tmp_Vz.nGrowVect(), WarpX::do_single_precision_comms, period);

    FillBoundary(tmp_Q_minus1_x, tmp_Q_minus1_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus2_x, tmp_Q_minus2_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus3_x, tmp_Q_minus3_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus4_x, tmp_Q_minus4_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus1_x, tmp_Q_plus1_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus2_x, tmp_Q_plus2_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus3_x, tmp_Q_plus3_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus4_x, tmp_Q_plus4_x.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus1_y, tmp_Q_minus1_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus2_y, tmp_Q_minus2_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus3_y, tmp_Q_minus3_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus4_y, tmp_Q_minus4_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus1_y, tmp_Q_plus1_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus2_y, tmp_Q_plus2_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus3_y, tmp_Q_plus3_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus4_y, tmp_Q_plus4_y.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus1_z, tmp_Q_minus1_z.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus2_z, tmp_Q_minus2_z.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus3_z, tmp_Q_minus3_z.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_minus4_z, tmp_Q_minus4_z.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus1_z, tmp_Q_plus1_z.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus2_z, tmp_Q_plus2_z.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus3_z, tmp_Q_plus3_z.nGrowVect(), WarpX::do_single_precision_comms, period);
    FillBoundary(tmp_Q_plus4_z, tmp_Q_plus4_z.nGrowVect(), WarpX::do_single_precision_comms, period);


        // Advection push
    #ifdef AMREX_USE_OMP
    #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    #endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box const &tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());
        amrex::Array4<Real> N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> NUz_arr = NU[lev][2]->array(mfi);

        amrex::Array4<amrex::Real> const &Vx = tmp_Vx.array(mfi);
        amrex::Array4<amrex::Real> const &Vy = tmp_Vy.array(mfi);
        amrex::Array4<amrex::Real> const &Vz = tmp_Vz.array(mfi);

        amrex::Array4<amrex::Real> const &Q_minus1_x = tmp_Q_minus1_x.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus2_x = tmp_Q_minus2_x.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus3_x = tmp_Q_minus3_x.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus4_x = tmp_Q_minus4_x.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus1_x = tmp_Q_plus1_x.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus2_x = tmp_Q_plus2_x.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus3_x = tmp_Q_plus3_x.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus4_x = tmp_Q_plus4_x.array(mfi);

        amrex::Array4<amrex::Real> const &Q_minus1_y = tmp_Q_minus1_y.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus2_y = tmp_Q_minus2_y.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus3_y = tmp_Q_minus3_y.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus4_y = tmp_Q_minus4_y.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus1_y = tmp_Q_plus1_y.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus2_y = tmp_Q_plus2_y.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus3_y = tmp_Q_plus3_y.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus4_y = tmp_Q_plus4_y.array(mfi);

        amrex::Array4<amrex::Real> const &Q_minus1_z = tmp_Q_minus1_z.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus2_z = tmp_Q_minus2_z.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus3_z = tmp_Q_minus3_z.array(mfi);
        amrex::Array4<amrex::Real> const &Q_minus4_z = tmp_Q_minus4_z.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus1_z = tmp_Q_plus1_z.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus2_z = tmp_Q_plus2_z.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus3_z = tmp_Q_plus3_z.array(mfi);
        amrex::Array4<amrex::Real> const &Q_plus4_z = tmp_Q_plus4_z.array(mfi);

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // Select the specific implmentation depending on dimensionality
                #if defined(WARPX_DIM_3D)


                // compute the fluxes:
                auto F1_minusx = flux(Q_minus1_x(i-1,j,k),Q_plus1_x(i,j,k),  Vx(i-1,j,k),Vx(i,j,k));
                auto F1_plusx =  flux(Q_minus1_x(i,j,k),  Q_plus1_x(i+1,j,k),Vx(i,j,k),  Vx(i+1,j,k));
                auto F2_minusx = flux(Q_minus2_x(i-1,j,k),Q_plus2_x(i,j,k),  Vx(i-1,j,k),Vx(i,j,k));
                auto F2_plusx =  flux(Q_minus2_x(i,j,k),  Q_plus2_x(i+1,j,k),Vx(i,j,k),  Vx(i+1,j,k));
                auto F3_minusx = flux(Q_minus3_x(i-1,j,k),Q_plus3_x(i,j,k),  Vx(i-1,j,k),Vx(i,j,k));
                auto F3_plusx =  flux(Q_minus3_x(i,j,k),  Q_plus3_x(i+1,j,k),Vx(i,j,k),  Vx(i+1,j,k));
                auto F4_minusx = flux(Q_minus4_x(i-1,j,k),Q_plus4_x(i,j,k),  Vx(i-1,j,k),Vx(i,j,k));
                auto F4_plusx =  flux(Q_minus4_x(i,j,k),  Q_plus4_x(i+1,j,k),Vx(i,j,k),  Vx(i+1,j,k));

                auto F1_minusy = flux(Q_minus1_y(i,j-1,k),Q_plus1_y(i,j,k),  Vy(i,j-1,k),Vy(i,j,k));
                auto F1_plusy =  flux(Q_minus1_y(i,j,k),  Q_plus1_y(i,j+1,k),Vy(i,j,k),  Vy(i,j+1,k));
                auto F2_minusy = flux(Q_minus2_y(i,j-1,k),Q_plus2_y(i,j,k),  Vy(i,j-1,k),Vy(i,j,k));
                auto F2_plusy =  flux(Q_minus2_y(i,j,k),  Q_plus2_y(i,j+1,k),Vy(i,j,k),  Vy(i,j+1,k));
                auto F3_minusy = flux(Q_minus3_y(i,j-1,k),Q_plus3_y(i,j,k),  Vy(i,j-1,k),Vy(i,j,k));
                auto F3_plusy =  flux(Q_minus3_y(i,j,k),  Q_plus3_y(i,j+1,k),Vy(i,j,k),  Vy(i,j+1,k));
                auto F4_minusy = flux(Q_minus4_y(i,j-1,k),Q_plus4_y(i,j,k),  Vy(i,j-1,k),Vy(i,j,k));
                auto F4_plusy =  flux(Q_minus4_y(i,j,k),  Q_plus4_y(i,j+1,k),Vy(i,j,k),  Vy(i,j+1,k));

                auto F1_minusz = flux(Q_minus1_z(i,j,k-1),Q_plus1_z(i,j,k),  Vz(i,j,k-1),Vz(i,j,k));
                auto F1_plusz =  flux(Q_minus1_z(i,j,k),  Q_plus1_z(i,j,k+1),Vz(i,j,k),  Vz(i,j,k+1));
                auto F2_minusz = flux(Q_minus2_z(i,j,k-1),Q_plus2_z(i,j,k),  Vz(i,j,k-1),Vz(i,j,k));
                auto F2_plusz =  flux(Q_minus2_z(i,j,k),  Q_plus2_z(i,j,k+1),Vz(i,j,k),  Vz(i,j,k+1));
                auto F3_minusz = flux(Q_minus3_z(i,j,k-1),Q_plus3_z(i,j,k),  Vz(i,j,k-1),Vz(i,j,k));
                auto F3_plusz =  flux(Q_minus3_z(i,j,k),  Q_plus3_z(i,j,k+1),Vz(i,j,k),  Vz(i,j,k+1));
                auto F4_minusz = flux(Q_minus4_z(i,j,k-1),Q_plus4_z(i,j,k),  Vz(i,j,k-1),Vz(i,j,k));
                auto F4_plusz =  flux(Q_minus4_z(i,j,k),  Q_plus4_z(i,j,k+1),Vz(i,j,k),  Vz(i,j,k+1));


                // Update Q from tn -> tn + dt
                N_arr(i,j,k) = N_arr(i,j,k) - cx*(F1_plusx - F1_minusx)
                                            - cy*(F1_plusy - F1_minusy)
                                            - cz*(F1_plusz - F1_minusz);
                NUx_arr(i,j,k) = NUx_arr(i,j,k) - cx*(F2_plusx - F2_minusx)
                                                - cy*(F2_plusy - F2_minusy)
                                                - cz*(F2_plusz - F2_minusz);
                NUy_arr(i,j,k) = NUy_arr(i,j,k) - cx*(F3_plusx - F3_minusx)
                                                - cy*(F3_plusy - F3_minusy)
                                                - cz*(F3_plusz - F3_minusz);
                NUz_arr(i,j,k) = NUz_arr(i,j,k) - cx*(F4_plusx - F4_minusx)
                                                - cy*(F4_plusy - F4_minusy)
                                                - cz*(F4_plusz - F4_minusz);

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

// Local helper functions for MUSCL-Handcock
// Rusanov Flux
amrex::Real WarpXFluidContainer::flux (amrex::Real Qm, amrex::Real Qp, amrex::Real Vm, amrex::Real Vp)
{
    auto c = std::max( std::abs(Vm) , std::abs(Vp) );
    return 0.5*(Vm*Qm + Vp*Qp) - (0.5*c)*(Qp - Qm);
}
// ave_minmod
amrex::Real WarpXFluidContainer::ave (amrex::Real a, amrex::Real b)
{
    if (a*b > 0.0)
        return minmod( maxmod(a,b), minmod(2.0*a,2.0*b));
    else
        return 0.0;
}
// mindmod
amrex::Real WarpXFluidContainer::minmod (amrex::Real a, amrex::Real b)
{
    if (a > 0.0 && b > 0.0)
        return std::min(a, b);
    else if (a < 0.0 && b < 0.0)
        return std::max(a, b);
    else
        return 0.0;
}
//maxmod
amrex::Real WarpXFluidContainer::maxmod (amrex::Real a, amrex::Real b)
{
    if (a > 0.0 && b > 0.0)
        return std::max(a, b);
    else if (a < 0.0 && b < 0.0)
        return std::min(a, b);
    else
        return 0.0;
}


// Momentum source from fields
void WarpXFluidContainer::GatherAndPush (
    int lev,
    const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
    const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz)
{

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
