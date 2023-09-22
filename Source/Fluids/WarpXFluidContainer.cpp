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
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXUtil.H"
using namespace ablastr::utils::communication;
using namespace amrex;

WarpXFluidContainer::WarpXFluidContainer(int nlevs_max, int ispecies, const std::string &name, const amrex::Geometry& geom)
{
    species_id = ispecies;
    species_name = name;

    plasma_injector = std::make_unique<PlasmaInjector>(species_id, species_name, geom);
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

        // default values of E_external and B_external
        // are used to set the E and B field when "constant" or "parser"
        // is not explicitly used in the input
        pp_species_name.query("B_ext_init_style", m_B_ext_s);
        std::transform(m_B_ext_s.begin(),
                       m_B_ext_s.end(),
                       m_B_ext_s.begin(),
                       ::tolower);
        pp_species_name.query("E_ext_init_style", m_E_ext_s);
        std::transform(m_E_ext_s.begin(),
                       m_E_ext_s.end(),
                       m_E_ext_s.begin(),
                       ::tolower);

        // Parse external fields:
        // if the input string for B_ext_s is
        // "parse_b_ext_function" then the mathematical expression
        // for the Bx_, By_, Bz_external_function(x,y,z)
        // must be provided in the input file.
        if (m_B_ext_s == "parse_b_ext_function") {
           // store the mathematical expression as string
           std::string str_Bx_ext_function;
           std::string str_By_ext_function;
           std::string str_Bz_ext_function;
           utils::parser::Store_parserString(
                pp_species_name, "Bx_external_function(x,y,z,t)",
                str_Bx_ext_function);
           utils::parser::Store_parserString(
                pp_species_name, "By_external_function(x,y,z,t)",
                str_By_ext_function);
           utils::parser::Store_parserString(
                pp_species_name, "Bz_external_function(x,y,z,t)",
                str_Bz_ext_function);

           // Parser for B_external on the fluid
           m_Bx_parser = std::make_unique<amrex::Parser>(
               utils::parser::makeParser(str_Bx_ext_function,{"x","y","z","t"}));
           m_By_parser = std::make_unique<amrex::Parser>(
               utils::parser::makeParser(str_By_ext_function,{"x","y","z","t"}));
           m_Bz_parser = std::make_unique<amrex::Parser>(
               utils::parser::makeParser(str_Bz_ext_function,{"x","y","z","t"}));

        }

        // if the input string for E_ext_s is
        // "parse_e_ext_function" then the mathematical expression
        // for the Ex_, Ey_, Ez_external_function(x,y,z)
        // must be provided in the input file.
        if (m_E_ext_s == "parse_e_ext_function") {
           // store the mathematical expression as string
           std::string str_Ex_ext_function;
           std::string str_Ey_ext_function;
           std::string str_Ez_ext_function;
           utils::parser::Store_parserString(
               pp_species_name, "Ex_external_function(x,y,z,t)",
               str_Ex_ext_function);
           utils::parser::Store_parserString(
               pp_species_name, "Ey_external_function(x,y,z,t)",
               str_Ey_ext_function);
           utils::parser::Store_parserString(
               pp_species_name, "Ez_external_function(x,y,z,t)",
               str_Ez_ext_function);
           // Parser for E_external on the fluid
           m_Ex_parser = std::make_unique<amrex::Parser>(
               utils::parser::makeParser(str_Ex_ext_function,{"x","y","z","t"}));
           m_Ey_parser = std::make_unique<amrex::Parser>(
               utils::parser::makeParser(str_Ey_ext_function,{"x","y","z","t"}));
           m_Ez_parser = std::make_unique<amrex::Parser>(
               utils::parser::makeParser(str_Ez_ext_function,{"x","y","z","t"}));

        }


        initialized = true;
    }
}

void WarpXFluidContainer::AllocateLevelMFs(int lev, const BoxArray &ba, const DistributionMapping &dm)
{
    int ncomps = 1;
    const amrex::IntVect nguards(AMREX_D_DECL(2, 2, 2));

    // set human-readable tag for each MultiFab
    auto const tag = [lev](std::string tagname)
    {
        tagname.append("[l=").append(std::to_string(lev)).append("]");
        return tagname;
    };

    WarpX::AllocInitMultiFab(N[lev], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid density"), 0.0_rt);

    WarpX::AllocInitMultiFab(NU[lev][0], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid momentum density [x]"), 0.0_rt);
    WarpX::AllocInitMultiFab(NU[lev][1], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid momentum density [y]"), 0.0_rt);
    WarpX::AllocInitMultiFab(NU[lev][2], amrex::convert(ba, amrex::IntVect::TheNodeVector()),
                            dm, ncomps, nguards, lev, tag("fluid momentum density [z]"), 0.0_rt);
}

void WarpXFluidContainer::InitData(int lev, amrex::Box init_box, amrex::Real cur_time)
{
    WARPX_PROFILE("WarpXFluidContainer::InitData");

    // Convert initialization box to nodal box
    init_box.surroundingNodes();

    // Extract objects that give the initial density and momentum
    InjectorDensity *inj_rho = plasma_injector->getInjectorDensity();
    InjectorMomentum *inj_mom = plasma_injector->getInjectorMomentumDevice();

    // Extract grid geometry properties
    WarpX &warpx = WarpX::GetInstance();
    const amrex::Geometry &geom = warpx.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const amrex::Real clight = PhysConst::c;
    const amrex::Real gamma_boost = WarpX::gamma_boost;
    const amrex::Real beta_boost = WarpX::beta_boost;

    // Loop through cells and initialize their value
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        amrex::Box tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());
        amrex::Array4<Real> const &N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> const &NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> const &NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> const &NUz_arr = NU[lev][2]->array(mfi);

        // Return the intersection of all cells and the ones we wish to update
        amrex::Box init_box_intersection = init_box & tile_box;

        amrex::ParallelFor(init_box_intersection,
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

                // Lorentz transform z (from boosted to lab frame)
                if (gamma_boost > 1._rt){
                    z = gamma_boost*(z + beta_boost*clight*cur_time);
                }

                amrex::Real n = inj_rho->getDensity(x, y, z);
                amrex::XDim3 u = inj_mom->getBulkMomentum(x, y, z);

                // Give u the right dimensions of m/s
                u.x = u.x * clight;
                u.y = u.y * clight;
                u.z = u.z * clight;

                // Check if n > 0 and if not, don't compute the boost
                // Lorentz transform n, u (from lab to boosted frame)
                if (n > 0.0){
                    if (gamma_boost > 1._rt){
                        amrex::Real gamma = sqrt(1.0 + (u.x*u.x + u.y*u.y + u.z*u.z)/(clight*clight));
                        amrex::Real n_boosted = gamma_boost*n*( 1.0 - beta_boost*u.z/(gamma*clight) );
                        amrex::Real uz_boosted = gamma_boost*(u.z - beta_boost*clight*gamma);
                        u.z = uz_boosted;
                        n = n_boosted;
                    }
                }

                // Multiply by clight so u is back in SI units
                N_arr(i, j, k) = n;
                NUx_arr(i, j, k) = n * u.x;
                NUy_arr(i, j, k) = n * u.y;
                NUz_arr(i, j, k) = n * u.z;

            }
        );
    }
}


void WarpXFluidContainer::Evolve(
    int lev,
    const amrex::MultiFab &Ex, const amrex::MultiFab &Ey, const amrex::MultiFab &Ez,
    const amrex::MultiFab &Bx, const amrex::MultiFab &By, const amrex::MultiFab &Bz,
    amrex::MultiFab* rho, amrex::MultiFab &jx, amrex::MultiFab &jy, amrex::MultiFab &jz,
    amrex::Real cur_time, bool skip_deposition)
{

    WARPX_PROFILE("WarpXFluidContainer::Evolve");

    if (rho && ! skip_deposition && ! do_not_deposit) {
         // Deposit charge before particle push, in component 0 of MultiFab rho.
         DepositCharge(lev, *rho, 0);
    }

    // Step the Lorentz Term
    if(!do_not_gather){
        GatherAndPush(lev, Ex, Ey, Ez, Bx, By, Bz, cur_time);
    }

    // Cylindrical centrifugal term
    if(!do_not_push){
#if defined(WARPX_DIM_RZ)
        centrifugal_source_rz(lev);
#endif

        // Apply (non-periodic) BC on the fluids (needed for spatial derivative),
        // and communicate N, NU at boundaries
        ApplyBcFluidsAndComms(lev);

        // Step the Advective term
        AdvectivePush_Muscl(lev);
    }

    // Deposit rho to the simulation mesh
    // Deposit charge (end of the step)
    if (rho && ! skip_deposition && ! do_not_deposit) {
        DepositCharge(lev, *rho, 1);
    }

    // Deposit J to the simulation mesh
    if (!skip_deposition && ! do_not_deposit) {
        DepositCurrent(lev, jx, jy, jz);
    }
}

// Momentum source due to curvature
void WarpXFluidContainer::ApplyBcFluidsAndComms (int lev)
{
    WARPX_PROFILE("WarpXFluidContainer::ApplyBcFluidsAndComms");

    WarpX &warpx = WarpX::GetInstance();
    const amrex::Geometry &geom = warpx.Geom(lev);
    const amrex::Periodicity &period = geom.periodicity();
    const Array<int,AMREX_SPACEDIM> periodic_directions = geom.isPeriodic();
    amrex::Box domain = geom.Domain();
    // Convert to nodal box
    domain.surroundingNodes();

    // H&C push the momentum
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        amrex::Box tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());

        amrex::Array4<Real> N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> NUz_arr = NU[lev][2]->array(mfi);

        //Grow the tilebox
        tile_box.grow(1);

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // If the cell is is first gaurd cell & the dimension is non
                // periodic, then copy Q_{i+1} = Q_{i-1}.
                // Don't check r-dir in Z:
#if defined(WARPX_DIM_3D)

                // Upper end (index 2)
                if ( (periodic_directions[2] != 1) && (k==domain.bigEnd(2)+1) ){
                    N_arr(i,j,k) = N_arr(i,j,k-2);
                    NUx_arr(i,j,k) = NUx_arr(i,j,k-2);
                    NUy_arr(i,j,k) = NUy_arr(i,j,k-2);
                    NUz_arr(i,j,k) = NUz_arr(i,j,k-2);

                // Lower end (index 2)
                } else if ( (periodic_directions[2] != 1) && (k==domain.smallEnd(2)-1) ) {
                    N_arr(i,j,k) = N_arr(i,j,k+2);
                    NUx_arr(i,j,k) = NUx_arr(i,j,k+2);
                    NUy_arr(i,j,k) = NUy_arr(i,j,k+2);
                    NUz_arr(i,j,k) = NUz_arr(i,j,k+2);
                }

#elif ( defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_3D) )

                // Upper end (index 1)
                if ( (periodic_directions[1] != 1) && (j==domain.bigEnd(1)+1) ){
                    N_arr(i,j,k) = N_arr(i,j-2,k);
                    NUx_arr(i,j,k) = NUx_arr(i,j-2,k);
                    NUy_arr(i,j,k) = NUy_arr(i,j-2,k);
                    NUz_arr(i,j,k) = NUz_arr(i,j-2,k);

                // Lower end (index 1`)
                } else if ( (periodic_directions[1] != 1) && (j==domain.smallEnd(1)-1) ) {
                    N_arr(i,j,k) = N_arr(i,j+2,k);
                    NUx_arr(i,j,k) = NUx_arr(i,j+2,k);
                    NUy_arr(i,j,k) = NUy_arr(i,j+2,k);
                    NUz_arr(i,j,k) = NUz_arr(i,j+2,k);

                }

#elif ( defined(WARPX_DIM_1D_Z) || defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D) )

                // Upper end (index 0)
                if ( (periodic_directions[0] != 1) && (i==domain.bigEnd(0)+1) ){
                    N_arr(i,j,k) = N_arr(i-2,j,k);
                    NUx_arr(i,j,k) = NUx_arr(i-2,j,k);
                    NUy_arr(i,j,k) = NUy_arr(i-2,j,k);
                    NUz_arr(i,j,k) = NUz_arr(i-2,j,k);

                // Lower end (index 0)
                } else if ( (periodic_directions[0] != 1) && (i==domain.smallEnd(0)-1) ) {
                    N_arr(i,j,k) = N_arr(i+2,j,k);
                    NUx_arr(i,j,k) = NUx_arr(i+2,j,k);
                    NUy_arr(i,j,k) = NUy_arr(i+2,j,k);
                    NUz_arr(i,j,k) = NUz_arr(i+2,j,k);
                }

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
#if defined(WARPX_DIM_3D)
    amrex::Real dt_over_dx = (dt/dx[0]);
    amrex::Real dt_over_dy = (dt/dx[1]);
    amrex::Real dt_over_dz = (dt/dx[2]);
    amrex::Real dt_over_dx_half = 0.5*(dt/dx[0]);
    amrex::Real dt_over_dy_half = 0.5*(dt/dx[1]);
    amrex::Real dt_over_dz_half = 0.5*(dt/dx[2]);
#elif defined(WARPX_DIM_XZ)
    amrex::Real dt_over_dx_half = 0.5*(dt/dx[0]);
    amrex::Real dt_over_dz_half = 0.5*(dt/dx[1]);
    amrex::Real dt_over_dx = (dt/dx[0]);
    amrex::Real dt_over_dz = (dt/dx[1]);
#elif defined(WARPX_DIM_RZ)
    const auto problo = geom.ProbLoArray();
    amrex::Real dt_over_dx_half = 0.5*(dt/dx[0]);
    amrex::Real dt_over_dz_half = 0.5*(dt/dx[1]);
    amrex::Box const& domain = geom.Domain();
#else
    amrex::Real dt_over_dz = (dt/dx[0]);
    amrex::Real dt_over_dz_half = 0.5*(dt/dx[0]);
#endif

    amrex::BoxArray ba = N[lev]->boxArray();

    // Temporary Half-step values
#if defined(WARPX_DIM_3D)
    amrex::MultiFab tmp_U_minus_x( amrex::convert(ba, IntVect(0,1,1)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_plus_x( amrex::convert(ba, IntVect(0,1,1)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_minus_y( amrex::convert(ba, IntVect(1,0,1)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_plus_y( amrex::convert(ba, IntVect(1,0,1)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_minus_z( amrex::convert(ba, IntVect(1,1,0)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_plus_z( amrex::convert(ba, IntVect(1,1,0)), N[lev]->DistributionMap(), 4, 1);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    amrex::MultiFab tmp_U_minus_x( amrex::convert(ba, IntVect(0,1)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_plus_x( amrex::convert(ba, IntVect(0,1)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_minus_z( amrex::convert(ba, IntVect(1,0)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_plus_z( amrex::convert(ba, IntVect(1,0)), N[lev]->DistributionMap(), 4, 1);
#else
    amrex::MultiFab tmp_U_minus_z( amrex::convert(ba, IntVect(0)), N[lev]->DistributionMap(), 4, 1);
    amrex::MultiFab tmp_U_plus_z( amrex::convert(ba, IntVect(0)), N[lev]->DistributionMap(), 4, 1);
#endif

    // Advection push
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        // Grow the entire domain
        amrex::Box box = mfi.validbox();
        box.grow(1);

        // Loop over a box with one extra gridpoint to avoid needing to communicate
        // the temporary arrays
        amrex::Box tile_box = mfi.growntilebox(1);

        // Limit the grown box for RZ at r = 0, r_max
#if defined (WARPX_DIM_RZ)
        const int idir = 0;
        const int n_cell = -1;
        tile_box.growLo(idir, n_cell);
        tile_box.growHi(idir, n_cell);
#endif

        amrex::Array4<Real> const &N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> const &NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> const &NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> const &NUz_arr = NU[lev][2]->array(mfi);

        // Only select tiles within the grown grid
#if defined(WARPX_DIM_3D)
        amrex::Box const box_x = amrex::convert( box, tmp_U_minus_x.ixType() );
        amrex::Box const box_y = amrex::convert( box, tmp_U_minus_y.ixType() );
        amrex::Box const box_z = amrex::convert( box, tmp_U_minus_z.ixType() );
        amrex::Array4<amrex::Real> U_minus_x = tmp_U_minus_x.array(mfi);
        amrex::Array4<amrex::Real> U_plus_x = tmp_U_plus_x.array(mfi);
        amrex::Array4<amrex::Real> U_minus_y = tmp_U_minus_y.array(mfi);
        amrex::Array4<amrex::Real> U_plus_y = tmp_U_plus_y.array(mfi);
        amrex::Array4<amrex::Real> U_minus_z = tmp_U_minus_z.array(mfi);
        amrex::Array4<amrex::Real> U_plus_z = tmp_U_plus_z.array(mfi);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        amrex::Box const box_x = amrex::convert( box, tmp_U_minus_x.ixType() );
        amrex::Box const box_z = amrex::convert( box, tmp_U_minus_z.ixType() );
        amrex::Array4<amrex::Real> U_minus_x = tmp_U_minus_x.array(mfi);
        amrex::Array4<amrex::Real> U_plus_x = tmp_U_plus_x.array(mfi);
        amrex::Array4<amrex::Real> U_minus_z = tmp_U_minus_z.array(mfi);
        amrex::Array4<amrex::Real> U_plus_z = tmp_U_plus_z.array(mfi);
#else
        amrex::Box const box_z = amrex::convert( box, tmp_U_minus_z.ixType() );
        amrex::Array4<amrex::Real> U_minus_z = tmp_U_minus_z.array(mfi);
        amrex::Array4<amrex::Real> U_plus_z = tmp_U_plus_z.array(mfi);
#endif

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                if( N_arr(i,j,k) > 0.0){

                    // - Grab local Uz Uy Ux gamma
                    // Isolate U from NU
                    amrex::Real Ux = (NUx_arr(i, j, k) / N_arr(i,j,k));
                    amrex::Real Uy = (NUy_arr(i, j, k) / N_arr(i,j,k));
                    amrex::Real Uz = (NUz_arr(i, j, k) / N_arr(i,j,k));

                    // Compute U_{x,y,z} at +/- 1 cell
#if defined(WARPX_DIM_3D)
                    amrex::Real Ux_px = 0.0, Ux_mx = 0.0, Uy_px = 0.0, Uy_mx = 0.0, Uz_px = 0.0, Uz_mx = 0.0;
                    if (N_arr(i+1,j,k) > 0.0) {
                        Ux_px = (NUx_arr(i+1, j, k) / N_arr(i+1,j,k));
                        Uy_px = (NUy_arr(i+1, j, k) / N_arr(i+1,j,k));
                        Uz_px = (NUz_arr(i+1, j, k) / N_arr(i+1,j,k));
                    }
                    if (N_arr(i-1,j,k) > 0.0) {
                        Ux_mx = (NUx_arr(i-1, j, k) / N_arr(i-1,j,k));
                        Uy_mx = (NUy_arr(i-1, j, k) / N_arr(i-1,j,k));
                        Uz_mx = (NUz_arr(i-1, j, k) / N_arr(i-1,j,k));
                    }
                    amrex::Real Ux_py = 0.0, Ux_my = 0.0, Uy_py = 0.0, Uy_my = 0.0, Uz_py = 0.0, Uz_my = 0.0;
                    if (N_arr(i,j+1,k) > 0.0) {
                        Ux_py = (NUx_arr(i, j+1, k) / N_arr(i,j+1,k));
                        Uy_py = (NUy_arr(i, j+1, k) / N_arr(i,j+1,k));
                        Uz_py = (NUz_arr(i, j+1, k) / N_arr(i,j+1,k));
                    }
                    if (N_arr(i,j-1,k) > 0.0) {
                        Ux_my = (NUx_arr(i, j-1, k) / N_arr(i,j-1,k));
                        Uy_my = (NUy_arr(i, j-1, k) / N_arr(i,j-1,k));
                        Uz_my = (NUz_arr(i, j-1, k) / N_arr(i,j-1,k));
                    }
                    amrex::Real Ux_pz = 0.0, Ux_mz = 0.0, Uy_pz = 0.0, Uy_mz = 0.0, Uz_pz = 0.0, Uz_mz = 0.0;
                    if (N_arr(i,j,k+1) > 0.0) {
                        Ux_pz = (NUx_arr(i, j, k+1) / N_arr(i,j,k+1));
                        Uy_pz = (NUy_arr(i, j, k+1) / N_arr(i,j,k+1));
                        Uz_pz = (NUz_arr(i, j, k+1) / N_arr(i,j,k+1));
                    }
                    if (N_arr(i,j,k-1) > 0.0) {
                        Ux_mz = (NUx_arr(i, j, k-1) / N_arr(i,j,k-1));
                        Uy_mz = (NUy_arr(i, j, k-1) / N_arr(i,j,k-1));
                        Uz_mz = (NUz_arr(i, j, k-1) / N_arr(i,j,k-1));
                    }

#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)

                    amrex::Real Ux_px = 0.0, Ux_mx = 0.0, Uy_px = 0.0, Uy_mx = 0.0, Uz_px = 0.0, Uz_mx = 0.0;
                    if (N_arr(i+1,j,k) > 0.0) {
                        Ux_px = (NUx_arr(i+1, j, k) / N_arr(i+1,j,k));
                        Uy_px = (NUy_arr(i+1, j, k) / N_arr(i+1,j,k));
                        Uz_px = (NUz_arr(i+1, j, k) / N_arr(i+1,j,k));
                    }
                    if (N_arr(i-1,j,k) > 0.0) {
                        Ux_mx = (NUx_arr(i-1, j, k) / N_arr(i-1,j,k));
                        Uy_mx = (NUy_arr(i-1, j, k) / N_arr(i-1,j,k));
                        Uz_mx = (NUz_arr(i-1, j, k) / N_arr(i-1,j,k));
                    }

                    amrex::Real Ux_pz = 0.0, Ux_mz = 0.0, Uy_pz = 0.0, Uy_mz = 0.0, Uz_pz = 0.0, Uz_mz = 0.0;
                    if (N_arr(i,j+1,k) > 0.0) {
                        Ux_pz = (NUx_arr(i, j+1, k) / N_arr(i,j+1,k));
                        Uy_pz = (NUy_arr(i, j+1, k) / N_arr(i,j+1,k));
                        Uz_pz = (NUz_arr(i, j+1, k) / N_arr(i,j+1,k));
                    }
                    if (N_arr(i,j-1,k) > 0.0) {
                        Ux_mz = (NUx_arr(i, j-1, k) / N_arr(i,j-1,k));
                        Uy_mz = (NUy_arr(i, j-1, k) / N_arr(i,j-1,k));
                        Uz_mz = (NUz_arr(i, j-1, k) / N_arr(i,j-1,k));
                    }

                // 1D
#else
                    amrex::Real Ux_pz = 0.0, Ux_mz = 0.0, Uy_pz = 0.0, Uy_mz = 0.0, Uz_pz = 0.0, Uz_mz = 0.0;
                    if (N_arr(i+1,j,k) > 0.0) {
                        Ux_pz = (NUx_arr(i+1, j, k) / N_arr(i+1,j,k));
                        Uy_pz = (NUy_arr(i+1, j, k) / N_arr(i+1,j,k));
                        Uz_pz = (NUz_arr(i+1, j, k) / N_arr(i+1,j,k));
                    }
                    if (N_arr(i-1,j,k) > 0.0) {
                        Ux_mz = (NUx_arr(i-1, j, k) / N_arr(i-1,j,k));
                        Uy_mz = (NUy_arr(i-1, j, k) / N_arr(i-1,j,k));
                        Uz_mz = (NUz_arr(i-1, j, k) / N_arr(i-1,j,k));
                    }
#endif


                    amrex::Real Uz_sq = Uz*Uz; amrex::Real Uy_sq = Uy*Uy; amrex::Real Ux_sq = Ux*Ux;
                    amrex::Real c_sq = clight*clight;
                    amrex::Real gamma = sqrt(1.0 + (Ux_sq + Uy_sq + Uz_sq)/(c_sq) );
                    amrex::Real a = c_sq*gamma*gamma*gamma;

                    // Calc Ax: (Needed for 2D, 3D, Rz)
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
                    amrex::Real Vx = Ux/gamma;
                    // Compute the Flux-Jacobian Elements in x
                    amrex::Real J00x = Vx;
                    amrex::Real J01x = N_arr(i,j,k)*(1/gamma)*(1-Vx*Vx/c_sq);
                    amrex::Real J02x = -N_arr(i,j,k)*Uy*Ux/a;
                    amrex::Real J03x = -N_arr(i,j,k)*Uz*Ux/a;

                    amrex::Real J10x = 0.0;
                    amrex::Real J11x = Vx;
                    amrex::Real J12x = 0.0;
                    amrex::Real J13x = 0.0;

                    amrex::Real J20x = 0.0;
                    amrex::Real J21x = 0.0;
                    amrex::Real J22x = Vx;
                    amrex::Real J23x = 0.0;

                    amrex::Real J30x = 0.0;
                    amrex::Real J31x = 0.0;
                    amrex::Real J32x = 0.0;
                    amrex::Real J33x = Vx;

#endif

                // Calc Ay: (Needed for 3d)
#if defined(WARPX_DIM_3D)
                    amrex::Real Vy = Uy/gamma;
                    // Compute the Flux-Jacobian Elements in y
                    amrex::Real J00y = Vy;
                    amrex::Real J01y = -N_arr(i,j,k)*Ux*Uy/a;
                    amrex::Real J02y = N_arr(i,j,k)*(1/gamma)*(1-Vy*Vy/c_sq);
                    amrex::Real J03y = -N_arr(i,j,k)*Uz*Uy/a;

                    amrex::Real J10y = 0.0;
                    amrex::Real J11y = Vy;
                    amrex::Real J12y = 0.0;
                    amrex::Real J13y = 0.0;

                    amrex::Real J20y = 0.0;
                    amrex::Real J21y = 0.0;
                    amrex::Real J22y = Vy;
                    amrex::Real J23y = 0.0;

                    amrex::Real J30y = 0.0;
                    amrex::Real J31y = 0.0;
                    amrex::Real J32y = 0.0;
                    amrex::Real J33y = Vy;


#endif
                    amrex::Real Vz = Uz/gamma;
                    // Calc Az: (needed for all)
                    // Compute the Flux-Jacobian Elements in z
                    amrex::Real J00z = Vz;
                    amrex::Real J01z = -N_arr(i,j,k)*Ux*Uz/a;
                    amrex::Real J02z = -N_arr(i,j,k)*Uy*Uz/a;
                    amrex::Real J03z = N_arr(i,j,k)*(1/gamma)*(1-Vz*Vz/c_sq);

                    amrex::Real J10z = 0.0;
                    amrex::Real J11z = Vz;
                    amrex::Real J12z = 0.0;
                    amrex::Real J13z = 0.0;

                    amrex::Real J20z = 0.0;
                    amrex::Real J21z = 0.0;
                    amrex::Real J22z = Vz;
                    amrex::Real J23z = 0.0;

                    amrex::Real J30z = 0.0;
                    amrex::Real J31z = 0.0;
                    amrex::Real J32z = 0.0;
                    amrex::Real J33z = Vz;


                // Select the specific implmentation depending on dimensionality
#if defined(WARPX_DIM_3D)

                    // Compute the cell slopes x
                    amrex::Real dU0x = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                    amrex::Real dU1x = ave( Ux - Ux_mx , Ux_px - Ux );
                    amrex::Real dU2x = ave( Uy - Uy_mx , Uy_px - Uy );
                    amrex::Real dU3x = ave( Uz - Uz_mx , Uz_px - Uz );

                    // Compute the cell slopes y
                    amrex::Real dU0y = ave( N_arr(i,j,k) - N_arr(i,j-1,k) , N_arr(i,j+1,k) - N_arr(i,j,k) );
                    amrex::Real dU1y = ave( Ux - Ux_my , Ux_py - Ux );
                    amrex::Real dU2y = ave( Uy - Uy_my , Uy_py - Uy );
                    amrex::Real dU3y = ave( Uz - Uz_my , Uz_py - Uz );

                    // Compute the cell slopes z
                    amrex::Real dU0z = ave( N_arr(i,j,k) - N_arr(i,j,k-1) , N_arr(i,j,k+1) - N_arr(i,j,k) );
                    amrex::Real dU1z = ave( Ux - Ux_mz , Ux_pz - Ux );
                    amrex::Real dU2z = ave( Uy - Uy_mz , Uy_pz - Uy );
                    amrex::Real dU3z = ave( Uz - Uz_mz , Uz_pz - Uz );

                    // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
                    amrex::Real JdU0x = J00x*dU0x + J01x*dU1x + J02x*dU2x + J03x*dU3x;
                    amrex::Real JdU1x = J10x*dU0x + J11x*dU1x + J12x*dU2x + J13x*dU3x;
                    amrex::Real JdU2x = J20x*dU0x + J21x*dU1x + J22x*dU2x + J23x*dU3x;
                    amrex::Real JdU3x = J30x*dU0x + J31x*dU1x + J32x*dU2x + J33x*dU3x;
                    amrex::Real JdU0y = J00y*dU0y + J01y*dU1y + J02y*dU2y + J03y*dU3y;
                    amrex::Real JdU1y = J10y*dU0y + J11y*dU1y + J12y*dU2y + J13y*dU3y;
                    amrex::Real JdU2y = J20y*dU0y + J21y*dU1y + J22y*dU2y + J23y*dU3y;
                    amrex::Real JdU3y = J30y*dU0y + J31y*dU1y + J32y*dU2y + J33y*dU3y;
                    amrex::Real JdU0z = J00z*dU0z + J01z*dU1z + J02z*dU2z + J03z*dU3z;
                    amrex::Real JdU1z = J10z*dU0z + J11z*dU1z + J12z*dU2z + J13z*dU3z;
                    amrex::Real JdU2z = J20z*dU0z + J21z*dU1z + J22z*dU2z + J23z*dU3z;
                    amrex::Real JdU3z = J30z*dU0z + J31z*dU1z + J32z*dU2z + J33z*dU3z;
                    amrex::Real U_tilde0 = N_arr(i,j,k)   - dt_over_dx_half*JdU0x - dt_over_dy_half*JdU0y - dt_over_dz_half*JdU0z;
                    amrex::Real U_tilde1 = Ux - dt_over_dx_half*JdU1x - dt_over_dy_half*JdU1y - dt_over_dz_half*JdU1z;
                    amrex::Real U_tilde2 = Uy - dt_over_dx_half*JdU2x - dt_over_dy_half*JdU2y - dt_over_dz_half*JdU2z;
                    amrex::Real U_tilde3 = Uz - dt_over_dx_half*JdU3x - dt_over_dy_half*JdU3y - dt_over_dz_half*JdU3z;


                    // Predict Q at the cell edges (x)
                    // (note that _plus is shifted due to grid location)
                    if ( box_x.contains(i,j,k) ) {
                        U_minus_x(i,j,k,0) = U_tilde0 + dU0x/2.0;
                        U_minus_x(i,j,k,1) = U_tilde1 + dU1x/2.0;
                        U_minus_x(i,j,k,2) = U_tilde2 + dU2x/2.0;
                        U_minus_x(i,j,k,3) = U_tilde3 + dU3x/2.0;
                    }
                    if ( box_x.contains(i-1,j,k) ) {
                        U_plus_x(i-1,j,k,0) = U_tilde0 - dU0x/2.0;
                        U_plus_x(i-1,j,k,1) = U_tilde1 - dU1x/2.0;
                        U_plus_x(i-1,j,k,2) = U_tilde2 - dU2x/2.0;
                        U_plus_x(i-1,j,k,3) = U_tilde3 - dU3x/2.0;
                    }

                    // Positivity and Monotonicty Limiter for density N:
                    if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) )) {
                        if ((U_minus_x(i,j,k,0) < 0.0) || (U_plus_x(i-1,j,k,0) < 0.0)){
                            U_minus_x(i,j,k,0) = N_arr(i,j,k);
                            U_minus_x(i,j,k,1) = Ux;
                            U_minus_x(i,j,k,2) = Uy;
                            U_minus_x(i,j,k,3) = Uz;
                            U_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                            U_plus_x(i-1,j,k,1) = Ux;
                            U_plus_x(i-1,j,k,2) = Uy;
                            U_plus_x(i-1,j,k,3) = Uz;
                        }
                    } else if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) != 1)) {
                        if (U_minus_x(i,j,k,0) < 0.0) {
                            U_minus_x(i,j,k,0) = N_arr(i,j,k);
                            U_minus_x(i,j,k,1) = Ux;
                            U_minus_x(i,j,k,2) = Uy;
                            U_minus_x(i,j,k,3) = Uz;
                        }
                    } else if (( box_x.contains(i,j,k) != 1 ) && ( box_x.contains(i-1,j,k) )) {
                        if (U_plus_x(i-1,j,k,0) < 0.0){
                            U_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                            U_plus_x(i-1,j,k,1) = Ux;
                            U_plus_x(i-1,j,k,2) = Uy;
                            U_plus_x(i-1,j,k,3) = Uz;
                        }
                    }

                    // Predict Q at the cell edges (y)
                    if ( box_y.contains(i,j,k) ) {
                        U_minus_y(i,j,k,0) = U_tilde0 + dU0y/2.0;
                        U_minus_y(i,j,k,1) = U_tilde1 + dU1y/2.0;
                        U_minus_y(i,j,k,2) = U_tilde2 + dU2y/2.0;
                        U_minus_y(i,j,k,3) = U_tilde3 + dU3y/2.0;
                    }
                    if ( box_y.contains(i,j-1,k) ) {
                        U_plus_y(i,j-1,k,0) = U_tilde0 - dU0y/2.0;
                        U_plus_y(i,j-1,k,1) = U_tilde1 - dU1y/2.0;
                        U_plus_y(i,j-1,k,2) = U_tilde2 - dU2y/2.0;
                        U_plus_y(i,j-1,k,3) = U_tilde3 - dU3y/2.0;
                    }

                    // Positivity and Monotonicty Limiter for density N:
                    if (( box_y.contains(i,j,k) ) && ( box_y.contains(i,j-1,k) )) {
                        if ((U_minus_y(i,j,k,0) < 0.0) || (U_plus_y(i,j-1,k,0) < 0.0)){
                            U_minus_y(i,j,k,0) = N_arr(i,j,k);
                            U_minus_y(i,j,k,1) = Ux;
                            U_minus_y(i,j,k,2) = Uy;
                            U_minus_y(i,j,k,3) = Uz;
                            U_plus_y(i,j-1,k,0) = N_arr(i,j,k);
                            U_plus_y(i,j-1,k,1) = Ux;
                            U_plus_y(i,j-1,k,2) = Uy;
                            U_plus_y(i,j-1,k,3) = Uz;
                        }
                    } else if (( box_y.contains(i,j,k) ) && ( box_y.contains(i,j-1,k) != 1)) {
                        if (U_minus_y(i,j,k,0) < 0.0) {
                            U_minus_y(i,j,k,0) = N_arr(i,j,k);
                            U_minus_y(i,j,k,1) = Ux;
                            U_minus_y(i,j,k,2) = Uy;
                            U_minus_y(i,j,k,3) = Uz;
                        }
                    } else if (( box_y.contains(i,j,k) != 1 ) && ( box_y.contains(i,j-1,k) )) {
                        if (U_plus_y(i,j-1,k,0) < 0.0){
                                U_plus_y(i,j-1,k,0) = N_arr(i,j,k);
                            U_plus_y(i,j-1,k,1) = Ux;
                            U_plus_y(i,j-1,k,2) = Uy;
                            U_plus_y(i,j-1,k,3) = Uz;
                        }
                    }

                    if ( box_z.contains(i,j,k) ) {
                    // Predict Q at the cell edges (z)
                        U_minus_z(i,j,k,0) = U_tilde0 + dU0z/2.0;
                        U_minus_z(i,j,k,1) = U_tilde1 + dU1z/2.0;
                        U_minus_z(i,j,k,2) = U_tilde2 + dU2z/2.0;
                        U_minus_z(i,j,k,3) = U_tilde3 + dU3z/2.0;
                    }
                    if ( box_z.contains(i,j,k-1) ) {
                        U_plus_z(i,j,k-1,0) = U_tilde0 - dU0z/2.0;
                        U_plus_z(i,j,k-1,1) = U_tilde1 - dU1z/2.0;
                        U_plus_z(i,j,k-1,2) = U_tilde2 - dU2z/2.0;
                        U_plus_z(i,j,k-1,3) = U_tilde3 - dU3z/2.0;
                    }


                    // Positivity and Monotonicty Limiter for density N: z
                    if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j,k-1) )) {
                        if ((U_minus_z(i,j,k,0) < 0.0) || (U_plus_z(i,j,k-1,0) < 0.0)){
                            U_minus_z(i,j,k,0) = N_arr(i,j,k);
                            U_minus_z(i,j,k,1) = Ux;
                            U_minus_z(i,j,k,2) = Uy;
                            U_minus_z(i,j,k,3) = Uz;
                            U_plus_z(i,j,k-1,0) = N_arr(i,j,k);
                            U_plus_z(i,j,k-1,1) = Ux;
                            U_plus_z(i,j,k-1,2) = Uy;
                            U_plus_z(i,j,k-1,3) = Uz;
                        }
                    } else if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j,k-1) != 1)) {
                        if (U_minus_z(i,j,k,0) < 0.0) {
                            U_minus_z(i,j,k,0) = N_arr(i,j,k);
                            U_minus_z(i,j,k,1) = Ux;
                            U_minus_z(i,j,k,2) = Uy;
                            U_minus_z(i,j,k,3) = Uz;
                        }
                    } else if (( box_z.contains(i,j,k) != 1 ) && ( box_z.contains(i,j,k-1) )) {
                        if (U_plus_z(i,j,k-1,0) < 0.0){
                            U_plus_z(i,j,k-1,0) = N_arr(i,j,k);
                            U_plus_z(i,j,k-1,1) = Ux;
                            U_plus_z(i,j,k-1,2) = Uy;
                            U_plus_z(i,j,k-1,3) = Uz;
                        }
                    }

#elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)

                    // Compute the cell slopes x
                    amrex::Real  dU0x = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                    amrex::Real  dU1x = ave( Ux - Ux_mx , Ux_px - Ux );
                    amrex::Real  dU2x = ave( Uy - Uy_mx , Uy_px - Uy );
                    amrex::Real  dU3x = ave( Uz - Uz_mx , Uz_px - Uz );

                    // Compute the cell slopes z
                    amrex::Real  dU0z = ave( N_arr(i,j,k) - N_arr(i,j-1,k) , N_arr(i,j+1,k) - N_arr(i,j,k) );
                    amrex::Real  dU1z = ave( Ux - Ux_mz , Ux_pz - Ux );
                    amrex::Real  dU2z = ave( Uy - Uy_mz , Uy_pz - Uy );
                    amrex::Real  dU3z = ave( Uz - Uz_mz , Uz_pz - Uz );
                    amrex::Real N_source = 0.0;

#if defined(WARPX_DIM_RZ)
                    amrex::Real dr = dx[0];
                    amrex::Real r = problo[0] + i * dr;
                    // Impose "none" boundaries
                    // Condition: dQx = 0 at r = 0
                    if  (i == domain.smallEnd(0)) {
                        // R|_{0+} -> L|_{0-}
                        // N -> N (N_arr(i-1,j,k) -> N_arr(i+1,j,k))
                        // NUr -> -NUr (NUx_arr(i-1,j,k) -> -NUx_arr(i+1,j,k))
                        // NUt -> -NUt (NUy_arr(i-1,j,k) -> -NUy_arr(i+1,j,k))
                        // NUz -> -NUz (NUz_arr(i-1,j,k) -> NUz_arr(i+1,j,k))
                        dU0x = ave( N_arr(i,j,k) - N_arr(i+1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                        dU1x = ave( Ux + Ux_px , Ux_px - Ux );
                        dU2x = ave( Uy + Uy_px , Uy_px - Uy );
                        dU3x = ave( Uz - Uz_px , Uz_px - Uz );
                    } else if (i == domain.bigEnd(0)+1) {
                        dU0z = ave( N_arr(i,j,k) - N_arr(i,j-1,k) , 0.0 );
                        dU1z = ave( Ux - Ux_mz , 0.0 );
                        dU2z = ave( Uy - Uy_mz , 0.0 );
                        dU3z = ave( Uz - Uz_mz , 0.0 );
                    }

                    // RZ sources:
                    if  (i != domain.smallEnd(0)) {
                        N_source = N_arr(i,j,k)*Vx/r;
                    }
#endif

                    // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
                    amrex::Real  JdU0x = J00x*dU0x + J01x*dU1x + J02x*dU2x + J03x*dU3x;
                    amrex::Real  JdU1x = J10x*dU0x + J11x*dU1x + J12x*dU2x + J13x*dU3x;
                    amrex::Real  JdU2x = J20x*dU0x + J21x*dU1x + J22x*dU2x + J23x*dU3x;
                    amrex::Real  JdU3x = J30x*dU0x + J31x*dU1x + J32x*dU2x + J33x*dU3x;
                    amrex::Real  JdU0z = J00z*dU0z + J01z*dU1z + J02z*dU2z + J03z*dU3z;
                    amrex::Real  JdU1z = J10z*dU0z + J11z*dU1z + J12z*dU2z + J13z*dU3z;
                    amrex::Real  JdU2z = J20z*dU0z + J21z*dU1z + J22z*dU2z + J23z*dU3z;
                    amrex::Real  JdU3z = J30z*dU0z + J31z*dU1z + J32z*dU2z + J33z*dU3z;
                    amrex::Real  U_tilde0 = N_arr(i,j,k)   - dt_over_dx_half*JdU0x - dt_over_dz_half*JdU0z - (dt/2.0)*N_source;
                    amrex::Real  U_tilde1 = Ux - dt_over_dx_half*JdU1x - dt_over_dz_half*JdU1z;
                    amrex::Real  U_tilde2 = Uy - dt_over_dx_half*JdU2x - dt_over_dz_half*JdU2z;
                    amrex::Real  U_tilde3 = Uz - dt_over_dx_half*JdU3x - dt_over_dz_half*JdU3z;

                    // Predict Q at the cell edges (x)
                    // (note that _plus is shifted due to grid location)
                    if ( box_x.contains(i,j,k) ) {
                        U_minus_x(i,j,k,0) = U_tilde0 + dU0x/2.0;
                        U_minus_x(i,j,k,1) = U_tilde1 + dU1x/2.0;
                        U_minus_x(i,j,k,2) = U_tilde2 + dU2x/2.0;
                        U_minus_x(i,j,k,3) = U_tilde3 + dU3x/2.0;
                    }
                    if ( box_x.contains(i-1,j,k) ) {
                        U_plus_x(i-1,j,k,0) = U_tilde0 - dU0x/2.0;
                        U_plus_x(i-1,j,k,1) = U_tilde1 - dU1x/2.0;
                        U_plus_x(i-1,j,k,2) = U_tilde2 - dU2x/2.0;
                        U_plus_x(i-1,j,k,3) = U_tilde3 - dU3x/2.0;
                    }

                // Positivity and Monotonicty Limiter for density N,
                // This sets the slope (dQ) to zero for all quantities
                    if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) )) {
                        if ((U_minus_x(i,j,k,0) < 0.0) || (U_plus_x(i-1,j,k,0) < 0.0)){
                            U_minus_x(i,j,k,0) = N_arr(i,j,k);
                            U_minus_x(i,j,k,1) = Ux;
                            U_minus_x(i,j,k,2) = Uy;
                            U_minus_x(i,j,k,3) = Uz;
                            U_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                            U_plus_x(i-1,j,k,1) = Ux;
                            U_plus_x(i-1,j,k,2) = Uy;
                            U_plus_x(i-1,j,k,3) = Uz;
                        }
                    } else if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) != 1)) {
                        if (U_minus_x(i,j,k,0) < 0.0) {
                            U_minus_x(i,j,k,0) = N_arr(i,j,k);
                            U_minus_x(i,j,k,1) = Ux;
                            U_minus_x(i,j,k,2) = Uy;
                            U_minus_x(i,j,k,3) = Uz;
                        }
                    } else if (( box_x.contains(i,j,k) != 1 ) && ( box_x.contains(i-1,j,k) )) {
                        if (U_plus_x(i-1,j,k,0) < 0.0){
                            U_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                            U_plus_x(i-1,j,k,1) = Ux;
                            U_plus_x(i-1,j,k,2) = Uy;
                            U_plus_x(i-1,j,k,3) = Uz;
                        }
                    }

                    if ( box_z.contains(i,j,k) ) {
                    // Predict Q at the cell edges (z)
                        U_minus_z(i,j,k,0) = U_tilde0 + dU0z/2.0;
                        U_minus_z(i,j,k,1) = U_tilde1 + dU1z/2.0;
                        U_minus_z(i,j,k,2) = U_tilde2 + dU2z/2.0;
                        U_minus_z(i,j,k,3) = U_tilde3 + dU3z/2.0;
                    }
                    if ( box_z.contains(i,j-1,k) ) {
                        U_plus_z(i,j-1,k,0) = U_tilde0 - dU0z/2.0;
                        U_plus_z(i,j-1,k,1) = U_tilde1 - dU1z/2.0;
                        U_plus_z(i,j-1,k,2) = U_tilde2 - dU2z/2.0;
                        U_plus_z(i,j-1,k,3) = U_tilde3 - dU3z/2.0;
                    }

                    // Positivity and Monotonicty Limiter for density N: z
                    if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j-1,k) )) {
                        if ((U_minus_z(i,j,k,0) < 0.0) || (U_plus_z(i,j-1,k,0) < 0.0)){
                            U_minus_z(i,j,k,0) = N_arr(i,j,k);
                            U_minus_z(i,j,k,1) = Ux;
                            U_minus_z(i,j,k,2) = Uy;
                            U_minus_z(i,j,k,3) = Uz;
                            U_plus_z(i,j-1,k,0) = N_arr(i,j,k);
                            U_plus_z(i,j-1,k,1) = Ux;
                            U_plus_z(i,j-1,k,2) = Uy;
                            U_plus_z(i,j-1,k,3) = Uz;
                        }
                    } else if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j-1,k) != 1)) {
                        if (U_minus_z(i,j,k,0) < 0.0) {
                            U_minus_z(i,j,k,0) = N_arr(i,j,k);
                            U_minus_z(i,j,k,1) = Ux;
                            U_minus_z(i,j,k,2) = Uy;
                            U_minus_z(i,j,k,3) = Uz;
                        }
                    } else if (( box_z.contains(i,j,k) != 1 ) && ( box_z.contains(i,j-1,k) )) {
                        if (U_plus_z(i,j-1,k,0) < 0.0){
                            U_plus_z(i,j-1,k,0) = N_arr(i,j,k);
                            U_plus_z(i,j-1,k,1) = Ux;
                            U_plus_z(i,j-1,k,2) = Uy;
                            U_plus_z(i,j-1,k,3) = Uz;
                        }
                    }

#else

                    // Compute the cell slopes z
                    amrex::Real  dU0z = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                    amrex::Real  dU1z = ave( Ux - Ux_mz , Ux_pz - Ux );
                    amrex::Real  dU2z = ave( Uy - Uy_mz , Uy_pz - Uy );
                    amrex::Real  dU3z = ave( Uz - Uz_mz , Uz_pz - Uz );

                    // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
                    amrex::Real  JdU0z = J00z*dU0z + J01z*dU1z + J02z*dU2z + J03z*dU3z;
                    amrex::Real  JdU1z = J10z*dU0z + J11z*dU1z + J12z*dU2z + J13z*dU3z;
                    amrex::Real  JdU2z = J20z*dU0z + J21z*dU1z + J22z*dU2z + J23z*dU3z;
                    amrex::Real  JdU3z = J30z*dU0z + J31z*dU1z + J32z*dU2z + J33z*dU3z;
                    amrex::Real  U_tilde0 = N_arr(i,j,k)   - dt_over_dz_half*JdU0z;
                    amrex::Real  U_tilde1 = Ux - dt_over_dz_half*JdU1z;
                    amrex::Real  U_tilde2 = Uy - dt_over_dz_half*JdU2z;
                    amrex::Real  U_tilde3 = Uz - dt_over_dz_half*JdU3z;

                    // Predict Q at the cell edges (z)
                    // (note that _plus is shifted due to grid location)
                    if ( box_z.contains(i,j,k) ) {
                    // Predict Q at the cell edges (z)
                        U_minus_z(i,j,k,0) = U_tilde0 + dU0z/2.0;
                        U_minus_z(i,j,k,1) = U_tilde1 + dU1z/2.0;
                        U_minus_z(i,j,k,2) = U_tilde2 + dU2z/2.0;
                        U_minus_z(i,j,k,3) = U_tilde3 + dU3z/2.0;
                    }
                    if ( box_z.contains(i-1,j,k) ) {
                        U_plus_z(i-1,j,k,0) = U_tilde0 - dU0z/2.0;
                        U_plus_z(i-1,j,k,1) = U_tilde1 - dU1z/2.0;
                        U_plus_z(i-1,j,k,2) = U_tilde2 - dU2z/2.0;
                        U_plus_z(i-1,j,k,3) = U_tilde3 - dU3z/2.0;
                    }

                    // Positivity and Monotonicty Limiter for density N,
                    // This sets the slope (dQ) to zero for all quantities
                    if (( box_z.contains(i,j,k) ) && ( box_z.contains(i-1,j,k) )) {
                        if ((U_minus_z(i,j,k,0) < 0.0) || (U_plus_z(i-1,j,k,0) < 0.0)) {
                            U_minus_z(i,j,k,0) = N_arr(i,j,k);
                            U_minus_z(i,j,k,1) = Ux;
                            U_minus_z(i,j,k,2) = Uy;
                            U_minus_z(i,j,k,3) = Uz;
                            U_plus_z(i-1,j,k,0) = N_arr(i,j,k);
                            U_plus_z(i-1,j,k,1) = Ux;
                            U_plus_z(i-1,j,k,2) = Uy;
                            U_plus_z(i-1,j,k,3) = Uz;
                        }
                    } else if (( box_z.contains(i,j,k) ) && ( box_z.contains(i-1,j,k) != 1)) {
                        if (U_minus_z(i,j,k,0) < 0.0) {
                            U_minus_z(i,j,k,0) = N_arr(i,j,k);
                            U_minus_z(i,j,k,1) = Ux;
                            U_minus_z(i,j,k,2) = Uy;
                            U_minus_z(i,j,k,3) = Uz;
                        }
                    } else if (( box_z.contains(i,j,k) != 1 ) && ( box_z.contains(i-1,j,k) )) {
                        if (U_plus_z(i-1,j,k,0) < 0.0){
                            U_plus_z(i-1,j,k,0) = N_arr(i,j,k);
                            U_plus_z(i-1,j,k,1) = Ux;
                            U_plus_z(i-1,j,k,2) = Uy;
                            U_plus_z(i-1,j,k,3) = Uz;
                        }
                    }

#endif
                // If N<= 0 then set the boundaries to zero
                } else {
#if defined(WARPX_DIM_3D) // 3D:
                    if ( box_x.contains(i,j,k) ) {
                        U_minus_x(i,j,k,0) = 0.0;
                        U_minus_x(i,j,k,1) = 0.0;
                        U_minus_x(i,j,k,2) = 0.0;
                        U_minus_x(i,j,k,3) = 0.0;
                    }
                    if ( box_x.contains(i-1,j,k) ) {
                        U_plus_x(i-1,j,k,0) = 0.0;
                        U_plus_x(i-1,j,k,1) = 0.0;
                        U_plus_x(i-1,j,k,2) = 0.0;
                        U_plus_x(i-1,j,k,3) = 0.0;
                    }
                    if ( box_y.contains(i,j,k) ) {
                        U_minus_y(i,j,k,0) = 0.0;
                        U_minus_y(i,j,k,1) = 0.0;
                        U_minus_y(i,j,k,2) = 0.0;
                        U_minus_y(i,j,k,3) = 0.0;
                    }
                    if ( box_y.contains(i,j-1,k) ) {
                        U_plus_y(i,j-1,k,0) = 0.0;
                        U_plus_y(i,j-1,k,1) = 0.0;
                        U_plus_y(i,j-1,k,2) = 0.0;
                        U_plus_y(i,j-1,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j,k) ) {
                        U_minus_z(i,j,k,0) = 0.0;
                        U_minus_z(i,j,k,1) = 0.0;
                        U_minus_z(i,j,k,2) = 0.0;
                        U_minus_z(i,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j,k-1) ) {
                        U_plus_z(i,j,k-1,0) = 0.0;
                        U_plus_z(i,j,k-1,1) = 0.0;
                        U_plus_z(i,j,k-1,2) = 0.0;
                        U_plus_z(i,j,k-1,3) = 0.0;
                    }
#elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ) // 2D:
                    if ( box_x.contains(i,j,k) ) {
                        U_minus_x(i,j,k,0) = 0.0;
                        U_minus_x(i,j,k,1) = 0.0;
                        U_minus_x(i,j,k,2) = 0.0;
                        U_minus_x(i,j,k,3) = 0.0;
                    }
                    if ( box_x.contains(i-1,j,k) ) {
                        U_plus_x(i-1,j,k,0) = 0.0;
                        U_plus_x(i-1,j,k,1) = 0.0;
                        U_plus_x(i-1,j,k,2) = 0.0;
                        U_plus_x(i-1,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j,k) ) {
                        U_minus_z(i,j,k,0) = 0.0;
                        U_minus_z(i,j,k,1) = 0.0;
                        U_minus_z(i,j,k,2) = 0.0;
                        U_minus_z(i,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j-1,k) ) {
                        U_plus_z(i,j-1,k,0) = 0.0;
                        U_plus_z(i,j-1,k,1) = 0.0;
                        U_plus_z(i,j-1,k,2) = 0.0;
                        U_plus_z(i,j-1,k,3) = 0.0;
                    }
#else // 1D:
                    if ( box_z.contains(i,j,k) ) {
                        U_minus_z(i,j,k,0) = 0.0;
                        U_minus_z(i,j,k,1) = 0.0;
                        U_minus_z(i,j,k,2) = 0.0;
                        U_minus_z(i,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i-1,j,k) ) {
                        U_plus_z(i-1,j,k,0) = 0.0;
                        U_plus_z(i-1,j,k,1) = 0.0;
                        U_plus_z(i-1,j,k,2) = 0.0;
                        U_plus_z(i-1,j,k,3) = 0.0;
                    }
#endif
                }
            }
        );
    }

    // Advection push
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box tile_box = mfi.tilebox(N[lev]->ixType().toIntVect());
        amrex::Array4<Real> N_arr = N[lev]->array(mfi);
        amrex::Array4<Real> NUx_arr = NU[lev][0]->array(mfi);
        amrex::Array4<Real> NUy_arr = NU[lev][1]->array(mfi);
        amrex::Array4<Real> NUz_arr = NU[lev][2]->array(mfi);

#if defined(WARPX_DIM_3D)
        amrex::Array4<amrex::Real> const &U_minus_x = tmp_U_minus_x.array(mfi);
        amrex::Array4<amrex::Real> const &U_plus_x = tmp_U_plus_x.array(mfi);
        amrex::Array4<amrex::Real> const &U_minus_y = tmp_U_minus_y.array(mfi);
        amrex::Array4<amrex::Real> const &U_plus_y = tmp_U_plus_y.array(mfi);
        amrex::Array4<amrex::Real> const &U_minus_z = tmp_U_minus_z.array(mfi);
        amrex::Array4<amrex::Real> const &U_plus_z = tmp_U_plus_z.array(mfi);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        amrex::Array4<amrex::Real> const &U_minus_x = tmp_U_minus_x.array(mfi);
        amrex::Array4<amrex::Real> const &U_plus_x = tmp_U_plus_x.array(mfi);
        amrex::Array4<amrex::Real> const &U_minus_z = tmp_U_minus_z.array(mfi);
        amrex::Array4<amrex::Real> const &U_plus_z = tmp_U_plus_z.array(mfi);
#else
        amrex::Array4<amrex::Real> const &U_minus_z = tmp_U_minus_z.array(mfi);
        amrex::Array4<amrex::Real> const &U_plus_z = tmp_U_plus_z.array(mfi);
#endif

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // Select the specific implmentation depending on dimensionality
#if defined(WARPX_DIM_3D)

                amrex::Real Vx_L_minus = 0.0, Vx_I_minus = 0.0, Vx_L_plus = 0.0, Vx_I_plus = 0.0;
                amrex::Real Vy_L_minus = 0.0, Vy_I_minus = 0.0, Vy_L_plus = 0.0, Vy_I_plus = 0.0;
                amrex::Real Vz_L_minus = 0.0, Vz_I_minus = 0.0, Vz_L_plus = 0.0, Vz_I_plus = 0.0;

                // Verify positive density, then compute velocity
                if (U_minus_x(i-1,j,k,0)>0.0) Vx_L_minus = V_calc(U_minus_x(i-1,j,k,1),U_minus_x(i-1,j,k,2),U_minus_x(i-1,j,k,3),clight,0);
                if (U_minus_x(i,j,k,0)>0.0)   Vx_I_minus = V_calc(U_minus_x(i,j,k,1),U_minus_x(i,j,k,2),U_minus_x(i,j,k,3),clight,0);
                if (U_plus_x(i-1,j,k,0)>0.0)   Vx_L_plus = V_calc(U_plus_x(i-1,j,k,1),U_plus_x(i-1,j,k,2),U_plus_x(i-1,j,k,3),clight,0);
                if (U_plus_x(i,j,k,0)>0.0) Vx_I_plus = V_calc(U_plus_x(i,j,k,1),U_plus_x(i,j,k,2),U_plus_x(i,j,k,3),clight,0);

                if (U_minus_y(i,j-1,k,0)>0.0) Vy_L_minus = V_calc(U_minus_y(i,j-1,k,1),U_minus_y(i,j-1,k,2),U_minus_y(i,j-1,k,3),clight,1);
                if (U_minus_y(i,j,k,0)>0.0)   Vy_I_minus = V_calc(U_minus_y(i,j,k,1),U_minus_y(i,j,k,2),U_minus_y(i,j,k,3),clight,1);
                if (U_plus_y(i,j-1,k,0)>0.0)   Vy_L_plus = V_calc(U_plus_y(i,j-1,k,1),U_plus_y(i,j-1,k,2),U_plus_y(i,j-1,k,3),clight,1);
                if (U_plus_y(i,j,k,0)>0.0)Vy_I_plus = V_calc(U_plus_y(i,j,k,1),U_plus_y(i,j,k,2),U_plus_y(i,j,k,3),clight,1);

                if (U_minus_z(i,j,k-1,0)>0.0) Vz_L_minus = V_calc(U_minus_z(i,j,k-1,1),U_minus_z(i,j,k-1,2),U_minus_z(i,j,k-1,3),clight,2);
                if (U_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(U_minus_z(i,j,k,1),U_minus_z(i,j,k,2),U_minus_z(i,j,k,3),clight,2);
                if (U_plus_z(i,j,k-1,0)>0.0)   Vz_L_plus = V_calc(U_plus_z(i,j,k-1,1),U_plus_z(i,j,k-1,2),U_plus_z(i,j,k-1,3),clight,2);
                if (U_plus_z(i,j,k,0)>0.0) Vz_I_plus = V_calc(U_plus_z(i,j,k,1),U_plus_z(i,j,k,2),U_plus_z(i,j,k,3),clight,2);

                // compute the fluxes:
                // (note that _plus is shifted due to grid location)
                amrex::Real F0_minusx = flux(U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F0_plusx =  flux(U_minus_x(i,j,k,0),  U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                amrex::Real F1_minusx = flux(U_minus_x(i-1,j,k,1)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,1)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F1_plusx =  flux(U_minus_x(i,j,k,1)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,1)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                amrex::Real F2_minusx = flux(U_minus_x(i-1,j,k,2)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,2)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F2_plusx =  flux(U_minus_x(i,j,k,2)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,2)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                amrex::Real F3_minusx = flux(U_minus_x(i-1,j,k,3)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,3)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F3_plusx =  flux(U_minus_x(i,j,k,3)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,3)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);

                amrex::Real F0_minusy = flux(U_minus_y(i,j-1,k,0),U_plus_y(i,j-1,k,0),  Vy_L_minus,Vy_L_plus);
                amrex::Real F0_plusy =  flux(U_minus_y(i,j,k,0),  U_plus_y(i,j,k,0),    Vy_I_minus,Vy_I_plus);
                amrex::Real F1_minusy = flux(U_minus_y(i,j-1,k,1)*U_minus_y(i,j-1,k,0),U_plus_y(i,j-1,k,1)*U_plus_y(i,j-1,k,0),  Vy_L_minus,Vy_L_plus);
                amrex::Real F1_plusy =  flux(U_minus_y(i,j,k,1)*U_minus_y(i,j,k,0),  U_plus_y(i,j,k,1)*U_plus_y(i,j,k,0),    Vy_I_minus,Vy_I_plus);
                amrex::Real F2_minusy = flux(U_minus_y(i,j-1,k,2)*U_minus_y(i,j-1,k,0),U_plus_y(i,j-1,k,2)*U_plus_y(i,j-1,k,0),  Vy_L_minus,Vy_L_plus);
                amrex::Real F2_plusy =  flux(U_minus_y(i,j,k,2)*U_minus_y(i,j,k,0),  U_plus_y(i,j,k,2)*U_plus_y(i,j,k,0),    Vy_I_minus,Vy_I_plus);
                amrex::Real F3_minusy = flux(U_minus_y(i,j-1,k,3)*U_minus_y(i,j-1,k,0),U_plus_y(i,j-1,k,3)*U_plus_y(i,j-1,k,0),  Vy_L_minus,Vy_L_plus);
                amrex::Real F3_plusy =  flux(U_minus_y(i,j,k,3)*U_minus_y(i,j,k,0),  U_plus_y(i,j,k,3)*U_plus_y(i,j,k,0),    Vy_I_minus,Vy_I_plus);

                amrex::Real F0_minusz = flux(U_minus_z(i,j,k-1,0),U_plus_z(i,j,k-1,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F0_plusz =  flux(U_minus_z(i,j,k,0),  U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F1_minusz = flux(U_minus_z(i,j,k-1,1)*U_minus_z(i,j,k-1,0),U_plus_z(i,j,k-1,1)*U_plus_z(i,j,k-1,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F1_plusz =  flux(U_minus_z(i,j,k,1)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,1)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F2_minusz = flux(U_minus_z(i,j,k-1,2)*U_minus_z(i,j,k-1,0),U_plus_z(i,j,k-1,2)*U_plus_z(i,j,k-1,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F2_plusz =  flux(U_minus_z(i,j,k,2)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,2)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F3_minusz = flux(U_minus_z(i,j,k-1,3)*U_minus_z(i,j,k-1,0),U_plus_z(i,j,k-1,3)*U_plus_z(i,j,k-1,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F3_plusz =  flux(U_minus_z(i,j,k,3)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,3)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);

                // Update Q from tn -> tn + dt
                N_arr(i,j,k) = N_arr(i,j,k) - dt_over_dx*(F0_plusx - F0_minusx)
                                            - dt_over_dy*(F0_plusy - F0_minusy)
                                            - dt_over_dz*(F0_plusz - F0_minusz);
                NUx_arr(i,j,k) = NUx_arr(i,j,k) - dt_over_dx*(F1_plusx - F1_minusx)
                                                - dt_over_dy*(F1_plusy - F1_minusy)
                                                - dt_over_dz*(F1_plusz - F1_minusz);
                NUy_arr(i,j,k) = NUy_arr(i,j,k) - dt_over_dx*(F2_plusx - F2_minusx)
                                                - dt_over_dy*(F2_plusy - F2_minusy)
                                                - dt_over_dz*(F2_plusz - F2_minusz);
                NUz_arr(i,j,k) = NUz_arr(i,j,k) - dt_over_dx*(F3_plusx - F3_minusx)
                                                - dt_over_dy*(F3_plusy - F3_minusy)
                                                - dt_over_dz*(F3_plusz - F3_minusz);

#elif defined(WARPX_DIM_XZ)

                amrex::Real Vx_L_minus = 0.0, Vx_I_minus = 0.0, Vx_L_plus = 0.0, Vx_I_plus = 0.0;
                amrex::Real Vz_L_minus = 0.0, Vz_I_minus = 0.0, Vz_L_plus = 0.0, Vz_I_plus = 0.0;

                // Verify positive density, then compute velocity
                if (U_minus_x(i-1,j,k,0)>0.0) Vx_L_minus = V_calc(U_minus_x(i-1,j,k,1),U_minus_x(i-1,j,k,2),U_minus_x(i-1,j,k,3),clight,0);
                if (U_minus_x(i,j,k,0)>0.0)   Vx_I_minus = V_calc(U_minus_x(i,j,k,1),U_minus_x(i,j,k,2),U_minus_x(i,j,k,3),clight,0);
                if (U_plus_x(i-1,j,k,0)>0.0)   Vx_L_plus = V_calc(U_plus_x(i-1,j,k,1),U_plus_x(i-1,j,k,2),U_plus_x(i-1,j,k,3),clight,0);
                if (U_plus_x(i,j,k,0)>0.0) Vx_I_plus = V_calc(U_plus_x(i,j,k,1),U_plus_x(i,j,k,2),U_plus_x(i,j,k,3),clight,0);

                if (U_minus_z(i,j-1,k,0)>0.0) Vz_L_minus = V_calc(U_minus_z(i,j-1,k,1),U_minus_z(i,j-1,k,2),U_minus_z(i,j-1,k,3),clight,2);
                if (U_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(U_minus_z(i,j,k,1),U_minus_z(i,j,k,2),U_minus_z(i,j,k,3),clight,2);
                if (U_plus_z(i,j-1,k,0)>0.0)   Vz_L_plus = V_calc(U_plus_z(i,j-1,k,1),U_plus_z(i,j-1,k,2),U_plus_z(i,j-1,k,3),clight,2);
                if (U_plus_z(i,j,k,0)>0.0) Vz_I_plus = V_calc(U_plus_z(i,j,k,1),U_plus_z(i,j,k,2),U_plus_z(i,j,k,3),clight,2);


                // compute the fluxes:
                // (note that _plus is shifted due to grid location)
                amrex::Real F0_minusx = flux(U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F0_plusx =  flux(U_minus_x(i,j,k,0),  U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                amrex::Real F1_minusx = flux(U_minus_x(i-1,j,k,1)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,1)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F1_plusx =  flux(U_minus_x(i,j,k,1)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,1)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                amrex::Real F2_minusx = flux(U_minus_x(i-1,j,k,2)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,2)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F2_plusx =  flux(U_minus_x(i,j,k,2)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,2)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                amrex::Real F3_minusx = flux(U_minus_x(i-1,j,k,3)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,3)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                amrex::Real F3_plusx =  flux(U_minus_x(i,j,k,3)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,3)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);

                amrex::Real F0_minusz = flux(U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F0_plusz =  flux(U_minus_z(i,j,k,0),  U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F1_minusz = flux(U_minus_z(i,j-1,k,1)*U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,1)*U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F1_plusz =  flux(U_minus_z(i,j,k,1)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,1)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F2_minusz = flux(U_minus_z(i,j-1,k,2)*U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,2)*U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F2_plusz =  flux(U_minus_z(i,j,k,2)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,2)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F3_minusz = flux(U_minus_z(i,j-1,k,3)*U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,3)*U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F3_plusz =  flux(U_minus_z(i,j,k,3)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,3)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);

                // Update Q from tn -> tn + dt
                N_arr(i,j,k) = N_arr(i,j,k) - dt_over_dx*(F0_plusx - F0_minusx)
                                            - dt_over_dz*(F0_plusz - F0_minusz);
                NUx_arr(i,j,k) = NUx_arr(i,j,k) - dt_over_dx*(F1_plusx - F1_minusx)
                                                - dt_over_dz*(F1_plusz - F1_minusz);
                NUy_arr(i,j,k) = NUy_arr(i,j,k) - dt_over_dx*(F2_plusx - F2_minusx)
                                                - dt_over_dz*(F2_plusz - F2_minusz);
                NUz_arr(i,j,k) = NUz_arr(i,j,k) - dt_over_dx*(F3_plusx - F3_minusx)
                                                - dt_over_dz*(F3_plusz - F3_minusz);

#elif defined(WARPX_DIM_RZ)

                // Compute the flux areas for RZ
                amrex::Real pi = 3.1415926535897932;
                // Cell-centered radius
                amrex::Real dr = dx[0];
                amrex::Real dz = dx[1]; // Must be 1
                amrex::Real r = problo[0] + i * dr;
                amrex::Real Vij = 0.0;
                amrex::Real S_Az = 0.0;

                // Volume element and z-facing surfaces
                if (i == domain.smallEnd(0)) {
                    Vij = 2.0*pi*(dr/2.0)*(dr/4.0)*dz;
                    S_Az = 2.0*pi*(dr/4.0)*(dr/2.0);
                } else if (i == domain.bigEnd(0)+1) {
                    Vij = 2.0*pi*(r - dr/4.0)*(dr/2.0)*dz;
                    S_Az = 2.0*pi*(r - dr/4.0)*(dr/2.0);
                }  else {
                    Vij = 2.0*pi*r*dr*dz;
                    S_Az = 2.0*pi*(r)*dr;
                }

                // Radial Surfaces
                amrex::Real S_Ar_plus = 2.0*pi*(r + dr/2.0)*dz;
                amrex::Real S_Ar_minus = 2.0*pi*(r - dr/2.0)*dz;
                if (i == domain.smallEnd(0))
                    S_Ar_minus = 0.0;
                if (i == domain.bigEnd(0)+1)
                    S_Ar_plus = 2.0*pi*(r)*dz;

                // Impose "none" boundaries
                // Condition: Vx(r) = 0 at boundaries
                amrex::Real Vx_L_minus = 0.0, Vx_I_minus = 0.0, Vx_L_plus = 0.0, Vx_I_plus = 0.0;
                amrex::Real Vz_L_minus = 0.0, Vz_I_minus = 0.0, Vz_L_plus = 0.0, Vz_I_plus = 0.0;
                if (U_minus_x(i,j,k,0)>0.0) Vx_I_minus = V_calc(U_minus_x(i,j,k,1),U_minus_x(i,j,k,2),U_minus_x(i,j,k,3),clight,0);
                if (U_plus_x(i-1,j,k,0)>0.0) Vx_L_plus = V_calc(U_plus_x(i-1,j,k,1),U_plus_x(i-1,j,k,2),U_plus_x(i-1,j,k,3),clight,0);

                if (U_minus_z(i,j-1,k,0)>0.0) Vz_L_minus = V_calc(U_minus_z(i,j-1,k,1),U_minus_z(i,j-1,k,2),U_minus_z(i,j-1,k,3),clight,2);
                if (U_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(U_minus_z(i,j,k,1),U_minus_z(i,j,k,2),U_minus_z(i,j,k,3),clight,2);
                if (U_plus_z(i,j-1,k,0)>0.0)   Vz_L_plus = V_calc(U_plus_z(i,j-1,k,1),U_plus_z(i,j-1,k,2),U_plus_z(i,j-1,k,3),clight,2);
                if (U_plus_z(i,j,k,0)>0.0) Vz_I_plus = V_calc(U_plus_z(i,j,k,1),U_plus_z(i,j,k,2),U_plus_z(i,j,k,3),clight,2);


                // compute the fluxes:
                // (note that _plus is shifted due to grid location)
                amrex::Real F0_minusx = 0.0, F1_minusx = 0.0, F2_minusx = 0.0, F3_minusx = 0.0;
                amrex::Real F0_plusx = 0.0, F1_plusx = 0.0, F2_plusx = 0.0, F3_plusx = 0.0;
                if (i != domain.smallEnd(0)) {
                    if (U_minus_x(i-1,j,k,0)>0.0) Vx_L_minus = V_calc(U_minus_x(i-1,j,k,1),U_minus_x(i-1,j,k,2),U_minus_x(i-1,j,k,3),clight,0);
                    F0_minusx = flux(U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                    F1_minusx = flux(U_minus_x(i-1,j,k,1)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,1)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                    F2_minusx = flux(U_minus_x(i-1,j,k,2)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,2)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                    F3_minusx = flux(U_minus_x(i-1,j,k,3)*U_minus_x(i-1,j,k,0),U_plus_x(i-1,j,k,3)*U_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                }
                if (i < domain.bigEnd(0)) {
                    if (U_plus_x(i,j,k,0)>0.0) Vx_I_plus = V_calc(U_plus_x(i,j,k,1),U_plus_x(i,j,k,2),U_plus_x(i,j,k,3),clight,0);
                    F0_plusx =  flux(U_minus_x(i,j,k,0),  U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                    F1_plusx =  flux(U_minus_x(i,j,k,1)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,1)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                    F2_plusx =  flux(U_minus_x(i,j,k,2)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,2)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                    F3_plusx =  flux(U_minus_x(i,j,k,3)*U_minus_x(i,j,k,0),  U_plus_x(i,j,k,3)*U_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                }

                amrex::Real F0_minusz = flux(U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus)*S_Az;
                amrex::Real F0_plusz =  flux(U_minus_z(i,j,k,0),  U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus)*S_Az;
                amrex::Real F1_minusz = flux(U_minus_z(i,j-1,k,1)*U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,1)*U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus)*S_Az;
                amrex::Real F1_plusz =  flux(U_minus_z(i,j,k,1)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,1)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus)*S_Az;
                amrex::Real F2_minusz = flux(U_minus_z(i,j-1,k,2)*U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,2)*U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus)*S_Az;
                amrex::Real F2_plusz =  flux(U_minus_z(i,j,k,2)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,2)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus)*S_Az;
                amrex::Real F3_minusz = flux(U_minus_z(i,j-1,k,3)*U_minus_z(i,j-1,k,0),U_plus_z(i,j-1,k,3)*U_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus)*S_Az;
                amrex::Real F3_plusz =  flux(U_minus_z(i,j,k,3)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,3)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus)*S_Az;

                // Update Q from tn -> tn + dt
                N_arr(i,j,k) = N_arr(i,j,k)     - (dt/Vij)*(F0_plusx - F0_minusx + F0_plusz - F0_minusz);
                NUx_arr(i,j,k) = NUx_arr(i,j,k) - (dt/Vij)*(F1_plusx - F1_minusx + F1_plusz - F1_minusz);
                NUy_arr(i,j,k) = NUy_arr(i,j,k) - (dt/Vij)*(F2_plusx - F2_minusx + F2_plusz - F2_minusz);
                NUz_arr(i,j,k) = NUz_arr(i,j,k) - (dt/Vij)*(F3_plusx - F3_minusx + F3_plusz - F3_minusz);

#else

                amrex::Real Vz_L_minus = 0.0, Vz_I_minus = 0.0, Vz_L_plus = 0.0, Vz_I_plus = 0.0;

                // Compute the half-timestep velocities
                if (U_minus_z(i-1,j,k,0)>0.0) Vz_L_minus = V_calc(U_minus_z(i-1,j,k,1),U_minus_z(i-1,j,k,2),U_minus_z(i-1,j,k,3),clight,2);
                if (U_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(U_minus_z(i,j,k,1),U_minus_z(i,j,k,2),U_minus_z(i,j,k,3),clight,2);
                if (U_plus_z(i-1,j,k,0)>0.0)  Vz_L_plus =  V_calc(U_plus_z(i-1,j,k,1),U_plus_z(i-1,j,k,2),U_plus_z(i-1,j,k,3),clight,2);
                if (U_plus_z(i,j,k,0)>0.0)    Vz_I_plus =  V_calc(U_plus_z(i,j,k,1),U_plus_z(i,j,k,2),U_plus_z(i,j,k,3),clight,2);

                // compute the fluxes:
                // (note that _plus is shifted due to grid location)
                amrex::Real F0_minusz = flux(U_minus_z(i-1,j,k,0),U_plus_z(i-1,j,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F0_plusz =  flux(U_minus_z(i,j,k,0),  U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F1_minusz = flux(U_minus_z(i-1,j,k,1)*U_minus_z(i-1,j,k,0),U_plus_z(i-1,j,k,1)*U_plus_z(i-1,j,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F1_plusz =  flux(U_minus_z(i,j,k,1)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,1)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F2_minusz = flux(U_minus_z(i-1,j,k,2)*U_minus_z(i-1,j,k,0),U_plus_z(i-1,j,k,2)*U_plus_z(i-1,j,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F2_plusz =  flux(U_minus_z(i,j,k,2)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,2)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                amrex::Real F3_minusz = flux(U_minus_z(i-1,j,k,3)*U_minus_z(i-1,j,k,0),U_plus_z(i-1,j,k,3)*U_plus_z(i-1,j,k,0),  Vz_L_minus,Vz_L_plus);
                amrex::Real F3_plusz =  flux(U_minus_z(i,j,k,3)*U_minus_z(i,j,k,0),  U_plus_z(i,j,k,3)*U_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);


                // Update Q from tn -> tn + dt
                N_arr(i,j,k) = N_arr(i,j,k)     - dt_over_dz*(F0_plusz - F0_minusz);
                NUx_arr(i,j,k) = NUx_arr(i,j,k) - dt_over_dz*(F1_plusz - F1_minusz);
                NUy_arr(i,j,k) = NUy_arr(i,j,k) - dt_over_dz*(F2_plusz - F2_minusz);
                NUz_arr(i,j,k) = NUz_arr(i,j,k) - dt_over_dz*(F3_plusz - F3_minusz);
#endif
            }
        );
    }
}


// Momentum source due to curvature
#if defined(WARPX_DIM_RZ)
void WarpXFluidContainer::centrifugal_source_rz (int lev)
{
    WARPX_PROFILE("WarpXFluidContainer::centrifugal_source_rz");

    WarpX &warpx = WarpX::GetInstance();
    const Real dt = warpx.getdt(lev);
    const amrex::Geometry &geom = warpx.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const amrex::Real clight = PhysConst::c;
    amrex::Box const& domain = geom.Domain();

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
        amrex::Array4<Real> const &NUz_arr = NU[lev][2]->array(mfi);

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // Verify density is non-zero
                if (N_arr(i,j,k)>0.0) {

                    // Compute r
                    amrex::Real r = problo[0] + i * dx[0];

                    // Isolate U from NU
                    amrex::Real u_r =     (NUx_arr(i, j, k) / (N_arr(i,j,k) * clight ));
                    amrex::Real u_theta = (NUy_arr(i, j, k) / (N_arr(i,j,k) * clight ));
                    amrex::Real u_z =     (NUz_arr(i, j, k) / (N_arr(i,j,k) * clight ));

                    // (SSP-RK3) Push the fluid momentum (R and Theta)
                    // F_r, F_theta are first order euler pushes of our rhs operator
                    if (i != domain.smallEnd(0)) {
                        amrex::Real u_r_1     = F_r(r,u_r,u_theta,u_z,dt);
                        amrex::Real u_theta_1 = F_theta(r,u_r,u_theta,u_z,dt);
                        amrex::Real u_r_2     = (0.75)*(u_r)     + (0.25)*F_r(r,u_r_1,u_theta_1,u_z,dt);
                        amrex::Real u_theta_2 = (0.75)*(u_theta) + (0.25)*F_theta(r,u_r_1,u_theta_1,u_z,dt);
                        u_r            = (1.0/3.0)*(u_r)     + (2.0/3.0)*F_r(r,u_r_2,u_theta_2,u_z,dt);
                        u_theta        = (1.0/3.0)*(u_theta) + (2.0/3.0)*F_theta(r,u_r_2,u_theta_2,u_z,dt);

                        // Calculate NU, save NUr, NUtheta
                        NUx_arr(i,j,k) = N_arr(i,j,k)*u_r*clight;
                        NUy_arr(i,j,k) = N_arr(i,j,k)*u_theta*clight;

                    // BC r = 0, u_theta = 0, and there is no extra source terms
                    } else {
                        NUx_arr(i,j,k) = 0.0;
                        NUy_arr(i,j,k) = 0.0;
                    }
                }
            }
        );
    }
}
#endif

// Momentum source from fields
void WarpXFluidContainer::GatherAndPush (
    int lev,
    const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
    const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
    Real t)
{
    WARPX_PROFILE("WarpXFluidContainer::GatherAndPush");

    WarpX &warpx = WarpX::GetInstance();
    const amrex::Real q = getCharge();
    const amrex::Real m = getMass();
    const Real dt = warpx.getdt(lev);
    const amrex::Geometry &geom = warpx.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const amrex::Real gamma_boost = WarpX::gamma_boost;
    const amrex::Real beta_boost = WarpX::beta_boost;
    //Check whether m_E_ext_s is "none"
    bool external_e_fields; // Needs intializing
    bool external_b_fields; // Needs intializing


   // Prepare interpolation of current components to cell center
    amrex::GpuArray<int, 3> Nodal_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> Ex_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> Ey_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> Ez_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> Bx_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> By_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> Bz_type = amrex::GpuArray<int, 3>{0, 0, 0};
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

    // External field parsers
    external_e_fields = (m_E_ext_s == "parse_e_ext_function");
    external_b_fields = (m_B_ext_s == "parse_b_ext_function");
    amrex::ParserExecutor<4> Exfield_parser;
    amrex::ParserExecutor<4> Eyfield_parser;
    amrex::ParserExecutor<4> Ezfield_parser;
    amrex::ParserExecutor<4> Bxfield_parser;
    amrex::ParserExecutor<4> Byfield_parser;
    amrex::ParserExecutor<4> Bzfield_parser;
    if (external_e_fields){
        constexpr int num_arguments = 4; //x,y,z,t
        Exfield_parser = m_Ex_parser->compile<num_arguments>();
        Eyfield_parser = m_Ey_parser->compile<num_arguments>();
        Ezfield_parser = m_Ez_parser->compile<num_arguments>();
    }

    if (external_b_fields){
        constexpr int num_arguments = 4; //x,y,z,t
        Bxfield_parser = m_Bx_parser->compile<num_arguments>();
        Byfield_parser = m_By_parser->compile<num_arguments>();
        Bzfield_parser = m_Bz_parser->compile<num_arguments>();
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
        amrex::GpuArray<int, 3U> coarsening_ratio = {1, 1, 1};

        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {

                // Only run if density is positive
                if (N_arr(i,j,k)>0.0) {

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

                    if (gamma_boost > 1._rt) { // Lorentz transform fields due to moving frame
                        if ( ( external_b_fields ) || ( external_e_fields ) ){

                            // Lorentz transform z (from boosted to lab frame)
                            amrex::Real Ex_ext_boost, Ey_ext_boost, Ez_ext_boost;
                            amrex::Real Bx_ext_boost, By_ext_boost, Bz_ext_boost;
                            amrex::Real Ex_ext_lab, Ey_ext_lab, Ez_ext_lab;
                            amrex::Real Bx_ext_lab, By_ext_lab, Bz_ext_lab;

                            // Grab the location
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

                            // Get the lab frame E and B
                            // Transform (boosted to lab)
                            amrex::Real t_lab = gamma_boost*(t + beta_boost*z/PhysConst::c);
                            amrex::Real z_lab = gamma_boost*(z + beta_boost*PhysConst::c*t);

                            // Grab the external fields in the lab frame:
                            if ( external_e_fields ) {
                                Ex_ext_lab = Exfield_parser(x, y, z_lab, t_lab);
                                Ey_ext_lab = Eyfield_parser(x, y, z_lab, t_lab);
                                Ez_ext_lab = Ezfield_parser(x, y, z_lab, t_lab);
                            }else{
                                Ex_ext_lab = 0.0;
                                Ey_ext_lab = 0.0;
                                Ez_ext_lab = 0.0;
                            }
                            if ( external_b_fields ) {
                                Bx_ext_lab = Bxfield_parser(x, y, z_lab, t_lab);
                                By_ext_lab = Byfield_parser(x, y, z_lab, t_lab);
                                Bz_ext_lab = Bzfield_parser(x, y, z_lab, t_lab);
                            }else{
                                Bx_ext_lab = 0.0;
                                By_ext_lab = 0.0;
                                Bz_ext_lab = 0.0;
                            }

                            // Transform E & B (lab to boosted frame)
                            // (Require both to for the lorentz transform)
                            // RHS m_parser
                            Ez_ext_boost = Ez_ext_lab;
                            Bz_ext_boost = Bz_ext_lab;
                            Ex_ext_boost = gamma_boost*(Ex_ext_lab - beta_boost*PhysConst::c*By_ext_lab);
                            Ey_ext_boost = gamma_boost*(Ey_ext_lab + beta_boost*PhysConst::c*Bx_ext_lab);
                            Bx_ext_boost = gamma_boost*(Bx_ext_lab + beta_boost*Ey_ext_lab/PhysConst::c);
                            By_ext_boost = gamma_boost*(By_ext_lab - beta_boost*Ex_ext_lab/PhysConst::c);

                            // Then add to Nodal quantities in the boosted frame:
                            Ex_Nodal += Ex_ext_boost;
                            Ey_Nodal += Ey_ext_boost;
                            Ez_Nodal += Ez_ext_boost;
                            Bx_Nodal += Bx_ext_boost;
                            By_Nodal += By_ext_boost;
                            Bz_Nodal += Bz_ext_boost;
                        }
                    } else {

                        // Added external e fields:
                        if ( external_e_fields ){
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

                            Ex_Nodal += Exfield_parser(x, y, z, t);
                            Ey_Nodal += Eyfield_parser(x, y, z, t);
                            Ez_Nodal += Ezfield_parser(x, y, z, t);
                        }

                        // Added external b fields:
                        if ( external_b_fields ){
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

                            Bx_Nodal += Bxfield_parser(x, y, z, t);
                            By_Nodal += Byfield_parser(x, y, z, t);
                            Bz_Nodal += Bzfield_parser(x, y, z, t);
                        }
                    }

                    // Isolate U from NU
                    amrex::Real tmp_Ux = (NUx_arr(i, j, k) / N_arr(i,j,k));
                    amrex::Real tmp_Uy = (NUy_arr(i, j, k) / N_arr(i,j,k));
                    amrex::Real tmp_Uz = (NUz_arr(i, j, k) / N_arr(i,j,k));

                    // Enforce RZ boundary conditions
#if defined(WARPX_DIM_RZ)
                    if  ( i == 0 ){
                        Ex_Nodal = 0.0;
                        Ey_Nodal = 0.0;
                        By_Nodal = 0.0;
                        Bx_Nodal = 0.0;
                    }
#endif

                    // Push the fluid momentum
                    UpdateMomentumHigueraCary(tmp_Ux, tmp_Uy, tmp_Uz,
                        Ex_Nodal, Ey_Nodal, Ez_Nodal,
                        Bx_Nodal, By_Nodal, Bz_Nodal, q, m, dt );

                    // Calculate NU
                    NUx_arr(i,j,k) = N_arr(i,j,k)*tmp_Ux;
                    NUy_arr(i,j,k) = N_arr(i,j,k)*tmp_Uy;
                    NUz_arr(i,j,k) = N_arr(i,j,k)*tmp_Uz;
                }
            }
        );
    }
}

void WarpXFluidContainer::DepositCharge (int lev, amrex::MultiFab &rho, int icomp)
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
        amrex::Array4<int> owner_mask_rho_arr = owner_mask_rho->array(mfi);

        // Deposit Rho
        amrex::ParallelFor(tile_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if ( owner_mask_rho_arr(i,j,k) ) rho_arr(i,j,k,icomp) += q*N_arr(i,j,k);
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
    amrex::GpuArray<int, 3> j_nodal_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> jx_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> jy_type = amrex::GpuArray<int, 3>{0, 0, 0};
    amrex::GpuArray<int, 3> jz_type = amrex::GpuArray<int, 3>{0, 0, 0};
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
                amrex::Real gamma = 1.0, Ux = 0.0, Uy = 0.0, Uz = 0.0;
                if (N_arr(i, j, k)>0.0){
                    Ux = NUx_arr(i, j, k)/N_arr(i, j, k);
                    Uy = NUy_arr(i, j, k)/N_arr(i, j, k);
                    Uz = NUz_arr(i, j, k)/N_arr(i, j, k);
                    gamma = std::sqrt(1.0 + ( Ux*Ux + Uy*Uy + Uz*Uz) * inv_clight_sq ) ;
                }
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

        amrex::Array4<int> owner_mask_x_arr = owner_mask_x->array(mfi);
        amrex::Array4<int> owner_mask_y_arr = owner_mask_y->array(mfi);
        amrex::Array4<int> owner_mask_z_arr = owner_mask_z->array(mfi);

        // When using the `Interp` function, one needs to specify whether coarsening is desired.
        // Here, we do not perform any coarsening.
        amrex::GpuArray<int, 3U> coarsening_ratio = {1, 1, 1};


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
