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

    // Loop through cells and initialize their value
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*N[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // TODO: Only run the code below if this cell overlaps is in `init_box`

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
                if (WarpX::gamma_boost > 1._rt){
                    z = WarpX::gamma_boost*(z + WarpX::beta_boost*clight*cur_time);
                }

                amrex::Real n = inj_rho->getDensity(x, y, z);
                auto u = inj_mom->getBulkMomentum(x, y, z);

                // Give u the right dimensions of m/s
                u.x = u.x * clight;
                u.y = u.y * clight;
                u.z = u.z * clight;

                // Check if n > 0 and if not, don't compute the boost
                // Lorentz transform n, u (from lab to boosted frame)
                if (n > 0.0){
                    if (WarpX::gamma_boost > 1._rt){
                        //amrex::Real n_boosted = WarpX::gamma_boost*(n - WarpX::beta_boost*n*u.z/clight);
                        //amrex::Real nuz_boosted = WarpX::gamma_boost*(n*u.z - WarpX::beta_boost*n*clight);
                        //u.z = nuz_boosted/n_boosted;
                        amrex::Real gamma = sqrt(1.0 + (u.x*u.x + u.y*u.y + u.z*u.z)/(clight*clight));
                        amrex::Real n_boosted = WarpX::gamma_boost*n*( 1.0 - WarpX::beta_boost*u.z/(gamma*clight) );
                        amrex::Real uz_boosted = WarpX::gamma_boost*(u.z - WarpX::beta_boost*clight*gamma);
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

    if (rho && ! skip_deposition && ! do_not_deposit) {
         // Deposit charge before particle push, in component 0 of MultiFab rho.
         DepositCharge(lev, *rho, 0);
    }

    // Step the Lorentz Term
    GatherAndPush(lev, Ex, Ey, Ez, Bx, By, Bz, cur_time);

    // Cylindrical centrifugal term
    #if defined(WARPX_DIM_RZ)
        centrifugal_source(lev);
    #endif

    // Apply (non-periodic) BC on the fluids (needed for spatial derivative),
    // and communicate N, NU at boundaries
    ApplyBcFluidsAndComms(lev);

    // Step the Advective term
    AdvectivePush_Muscl(lev);

    // Deposit rho to the simulation mesh
    // Deposit charge (end of the step)
    if (rho && ! skip_deposition && ! do_not_deposit) {
        DepositCharge(lev, *rho, 1);
    }

    // Deposit J to the simulation mesh
    if (!skip_deposition) {
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
        auto cx = (dt/dx[0]);
        auto cy = (dt/dx[1]);
        auto cz = (dt/dx[2]);
        auto cx_half = 0.5*(dt/dx[0]);
        auto cy_half = 0.5*(dt/dx[1]);
        auto cz_half = 0.5*(dt/dx[2]);
    #elif defined(WARPX_DIM_XZ)
        auto cx_half = 0.5*(dt/dx[0]);
        auto cz_half = 0.5*(dt/dx[1]);
        auto cx = (dt/dx[0]);
        auto cz = (dt/dx[1]);
    #elif defined(WARPX_DIM_RZ)
        const auto problo = geom.ProbLoArray();
        auto cx_half = 0.5*(dt/dx[0]);
        auto cz_half = 0.5*(dt/dx[1]);
        amrex::Box const& domain = geom.Domain();
    #else
        auto cz = (dt/dx[0]);
        auto cz_half = 0.5*(dt/dx[0]);
    #endif

    amrex::BoxArray ba = N[lev]->boxArray();

    // Temporary Half-step values
    #if defined(WARPX_DIM_3D)
        amrex::MultiFab tmp_Q_minus_x( amrex::convert(ba, IntVect(0,1,1)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_plus_x( amrex::convert(ba, IntVect(0,1,1)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_minus_y( amrex::convert(ba, IntVect(1,0,1)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_plus_y( amrex::convert(ba, IntVect(1,0,1)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_minus_z( amrex::convert(ba, IntVect(1,1,0)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_plus_z( amrex::convert(ba, IntVect(1,1,0)), N[lev]->DistributionMap(), 4, 1);
    #elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        amrex::MultiFab tmp_Q_minus_x( amrex::convert(ba, IntVect(0,1)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_plus_x( amrex::convert(ba, IntVect(0,1)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_minus_z( amrex::convert(ba, IntVect(1,0)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_plus_z( amrex::convert(ba, IntVect(1,0)), N[lev]->DistributionMap(), 4, 1);
    #else
        amrex::MultiFab tmp_Q_minus_z( amrex::convert(ba, IntVect(0)), N[lev]->DistributionMap(), 4, 1);
        amrex::MultiFab tmp_Q_plus_z( amrex::convert(ba, IntVect(0)), N[lev]->DistributionMap(), 4, 1);
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
            amrex::Box const box_x = amrex::convert( box, tmp_Q_minus_x.ixType() );
            amrex::Box const box_y = amrex::convert( box, tmp_Q_minus_y.ixType() );
            amrex::Box const box_z = amrex::convert( box, tmp_Q_minus_z.ixType() );
            amrex::Array4<amrex::Real> Q_minus_x = tmp_Q_minus_x.array(mfi);
            amrex::Array4<amrex::Real> Q_plus_x = tmp_Q_plus_x.array(mfi);
            amrex::Array4<amrex::Real> Q_minus_y = tmp_Q_minus_y.array(mfi);
            amrex::Array4<amrex::Real> Q_plus_y = tmp_Q_plus_y.array(mfi);
            amrex::Array4<amrex::Real> Q_minus_z = tmp_Q_minus_z.array(mfi);
            amrex::Array4<amrex::Real> Q_plus_z = tmp_Q_plus_z.array(mfi);
        #elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            amrex::Box const box_x = amrex::convert( box, tmp_Q_minus_x.ixType() );
            amrex::Box const box_z = amrex::convert( box, tmp_Q_minus_z.ixType() );
            amrex::Array4<amrex::Real> Q_minus_x = tmp_Q_minus_x.array(mfi);
            amrex::Array4<amrex::Real> Q_plus_x = tmp_Q_plus_x.array(mfi);
            amrex::Array4<amrex::Real> Q_minus_z = tmp_Q_minus_z.array(mfi);
            amrex::Array4<amrex::Real> Q_plus_z = tmp_Q_plus_z.array(mfi);
        #else
            amrex::Box const box_z = amrex::convert( box, tmp_Q_minus_z.ixType() );
            amrex::Array4<amrex::Real> Q_minus_z = tmp_Q_minus_z.array(mfi);
            amrex::Array4<amrex::Real> Q_plus_z = tmp_Q_plus_z.array(mfi);
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
                    amrex::Real Uz_sq = Uz*Uz; amrex::Real Uy_sq = Uy*Uy; amrex::Real Ux_sq = Ux*Ux;
                    amrex::Real Uz_cubed = Uz_sq*Uz; amrex::Real Uy_cubed = Uy_sq*Uy; amrex::Real Ux_cubed = Ux_sq*Ux;
                    amrex::Real c_sq = clight*clight;
                    amrex::Real gamma = sqrt(1.0 + (Ux_sq + Uy_sq + Uz_sq)/(c_sq) );
                    amrex::Real gamma_cubed = gamma*gamma*gamma;
                    amrex::Real a = c_sq*gamma_cubed;

                    // Calc Ax: (Needed for 2D, 3D, Rz)
                    #if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
                        // Compute the Flux-Jacobian Elements in x
                        amrex::Real A00x = (Ux*(Uz_sq)+Ux*(Uy_sq)+(Ux_cubed))/a;
                        amrex::Real A01x = ((c_sq)+(Uz_sq)+(Uy_sq))/a;
                        amrex::Real A02x = -(Ux*Uy)/a;
                        amrex::Real A03x = -(Ux*Uz)/a;

                        amrex::Real A10x = -(Ux_sq)/(gamma_cubed);
                        amrex::Real A11x = (2.0*Ux*(c_sq)+2.0*Ux*(Uz_sq)+2.0*Ux*(Uy_sq)+(Ux_cubed))/a;
                        amrex::Real A12x = -((Ux_sq)*Uy)/a;
                        amrex::Real A13x = -((Ux_sq)*Uz)/a;

                        amrex::Real A20x = -(Ux*Uy)/(gamma_cubed);
                        amrex::Real A21x = (Uy*(c_sq)+Uy*(Uz_sq)+(Uy_cubed))/a;
                        amrex::Real A22x = (Ux*(c_sq)+Ux*(Uz_sq)+(Ux_cubed))/a;
                        amrex::Real A23x = -(Ux*Uy*Uz)/a;

                        amrex::Real A30x = -(Ux*Uz)/(gamma_cubed);
                        amrex::Real A31x = (Uz*(c_sq)+(Uz_cubed)+(Uy_sq)*Uz)/a;
                        amrex::Real A32x = -(Ux*Uy*Uz)/a;
                        amrex::Real A33x = (Ux*(c_sq)+Ux*(Uy_sq)+(Ux_cubed))/a;
                    #endif

                    // Calc Ay: (Needed for 3d)
                    #if defined(WARPX_DIM_3D)
                        // Compute the Flux-Jacobian Elements in y
                        amrex::Real A00y = (Uy*(Uz_sq)+(Uy_cubed)+(Ux_sq)*Uy)/a;
                        amrex::Real A01y = -(Ux*Uy)/a;
                        amrex::Real A02y = ((c_sq)+(Uz_sq)+(Ux_sq))/a;
                        amrex::Real A03y = -(Uy*Uz)/a;

                        amrex::Real A10y = -(Ux*Uy)/(gamma_cubed);
                        amrex::Real A11y = (Uy*(c_sq)+Uy*(Uz_sq)+(Uy_cubed))/a;
                        amrex::Real A12y = (Ux*(c_sq)+Ux*(Uz_sq)+(Ux_cubed))/a;
                        amrex::Real A13y = -(Ux*Uy*Uz)/a;

                        amrex::Real A20y = -(Uy_sq)/(gamma_cubed);
                        amrex::Real A21y = -(Ux*(Uy_sq))/a;
                        amrex::Real A22y = (2.0*Uy*(c_sq)+2.0*Uy*(Uz_sq)+(Uy_cubed)+2.0*(Ux_sq)*Uy)/a;
                        amrex::Real A23y = -((Uy_sq)*Uz)/a;

                        amrex::Real A30y = -(Uy*Uz)/(gamma_cubed);
                        amrex::Real A31y = -(Ux*Uy*Uz)/a;
                        amrex::Real A32y = (Uz*(c_sq)+(Uz_cubed)+(Ux_sq)*Uz)/a;
                        amrex::Real A33y = (Uy*(c_sq)+(Uy_cubed)+(Ux_sq)*Uy)/a;

                    #endif

                        // Calc Az: (needed for all)
                        // Compute the Flux-Jacobian Elements in z
                        amrex::Real A00z = ((Uz_cubed)+((Uy_sq)+(Ux_sq))*Uz)/a;
                        amrex::Real A01z = -(Ux*Uz)/a;
                        amrex::Real A02z = -(Uy*Uz)/a;
                        amrex::Real A03z = ((c_sq)+(Uy_sq)+(Ux_sq))/a;

                        amrex::Real A10z = -(Ux*Uz)/(gamma_cubed);
                        amrex::Real A11z = (Uz*(c_sq)+(Uz_cubed)+(Uy_sq)*Uz)/a;
                        amrex::Real A12z = -(Ux*Uy*Uz)/a;
                        amrex::Real A13z = (Ux*(c_sq)+Ux*(Uy_sq)+(Ux_cubed))/a;

                        amrex::Real A20z = -(Uy*Uz)/(gamma_cubed);
                        amrex::Real A21z = -(Ux*Uy*Uz)/a;
                        amrex::Real A22z = (Uz*(c_sq)+(Uz_cubed)+(Ux_sq)*Uz)/a;
                        amrex::Real A23z = (Uy*(c_sq)+(Uy_cubed)+(Ux_sq)*Uy)/a;

                        amrex::Real A30z = -(Uz_sq)/(gamma_cubed);
                        amrex::Real A31z = -(Ux*(Uz_sq))/a;
                        amrex::Real A32z = -(Uy*(Uz_sq))/a;
                        amrex::Real A33z = (2.0*Uz*(c_sq)+(Uz_cubed)+(2.0*(Uy_sq)+2.0*(Ux_sq))*Uz)/a;


                    // Select the specific implmentation depending on dimensionality
                    #if defined(WARPX_DIM_3D)

                        // Compute the cell slopes x
                        amrex::Real dQ0x = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                        amrex::Real dQ1x = ave( NUx_arr(i,j,k) - NUx_arr(i-1,j,k) , NUx_arr(i+1,j,k) - NUx_arr(i,j,k) );
                        amrex::Real dQ2x = ave( NUy_arr(i,j,k) - NUy_arr(i-1,j,k) , NUy_arr(i+1,j,k) - NUy_arr(i,j,k) );
                        amrex::Real dQ3x = ave( NUz_arr(i,j,k) - NUz_arr(i-1,j,k) , NUz_arr(i+1,j,k) - NUz_arr(i,j,k) );

                        // Compute the cell slopes y
                        amrex::Real dQ0y = ave( N_arr(i,j,k) - N_arr(i,j-1,k) , N_arr(i,j+1,k) - N_arr(i,j,k) );
                        amrex::Real dQ1y = ave( NUx_arr(i,j,k) - NUx_arr(i,j-1,k) , NUx_arr(i,j+1,k) - NUx_arr(i,j,k) );
                        amrex::Real dQ2y = ave( NUy_arr(i,j,k) - NUy_arr(i,j-1,k) , NUy_arr(i,j+1,k) - NUy_arr(i,j,k) );
                        amrex::Real dQ3y = ave( NUz_arr(i,j,k) - NUz_arr(i,j-1,k) , NUz_arr(i,j+1,k) - NUz_arr(i,j,k) );

                        // Compute the cell slopes z
                        amrex::Real dQ0z = ave( N_arr(i,j,k) - N_arr(i,j,k-1) , N_arr(i,j,k+1) - N_arr(i,j,k) );
                        amrex::Real dQ1z = ave( NUx_arr(i,j,k) - NUx_arr(i,j,k-1) , NUx_arr(i,j,k+1) - NUx_arr(i,j,k) );
                        amrex::Real dQ2z = ave( NUy_arr(i,j,k) - NUy_arr(i,j,k-1) , NUy_arr(i,j,k+1) - NUy_arr(i,j,k) );
                        amrex::Real dQ3z = ave( NUz_arr(i,j,k) - NUz_arr(i,j,k-1) , NUz_arr(i,j,k+1) - NUz_arr(i,j,k) );

                        // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
                        amrex::Real AdQ0x = A00x*dQ0x + A01x*dQ1x + A02x*dQ2x + A03x*dQ3x;
                        amrex::Real AdQ1x = A10x*dQ0x + A11x*dQ1x + A12x*dQ2x + A13x*dQ3x;
                        amrex::Real AdQ2x = A20x*dQ0x + A21x*dQ1x + A22x*dQ2x + A23x*dQ3x;
                        amrex::Real AdQ3x = A30x*dQ0x + A31x*dQ1x + A32x*dQ2x + A33x*dQ3x;
                        amrex::Real AdQ0y = A00y*dQ0y + A01y*dQ1y + A02y*dQ2y + A03y*dQ3y;
                        amrex::Real AdQ1y = A10y*dQ0y + A11y*dQ1y + A12y*dQ2y + A13y*dQ3y;
                        amrex::Real AdQ2y = A20y*dQ0y + A21y*dQ1y + A22y*dQ2y + A23y*dQ3y;
                        amrex::Real AdQ3y = A30y*dQ0y + A31y*dQ1y + A32y*dQ2y + A33y*dQ3y;
                        amrex::Real AdQ0z = A00z*dQ0z + A01z*dQ1z + A02z*dQ2z + A03z*dQ3z;
                        amrex::Real AdQ1z = A10z*dQ0z + A11z*dQ1z + A12z*dQ2z + A13z*dQ3z;
                        amrex::Real AdQ2z = A20z*dQ0z + A21z*dQ1z + A22z*dQ2z + A23z*dQ3z;
                        amrex::Real AdQ3z = A30z*dQ0z + A31z*dQ1z + A32z*dQ2z + A33z*dQ3z;
                        amrex::Real Q_tilde0 = N_arr(i,j,k)   - cx_half*AdQ0x - cy_half*AdQ0y - cz_half*AdQ0z;
                        amrex::Real Q_tilde1 = NUx_arr(i,j,k) - cx_half*AdQ1x - cy_half*AdQ1y - cz_half*AdQ1z;
                        amrex::Real Q_tilde2 = NUy_arr(i,j,k) - cx_half*AdQ2x - cy_half*AdQ2y - cz_half*AdQ2z;
                        amrex::Real Q_tilde3 = NUz_arr(i,j,k) - cx_half*AdQ3x - cy_half*AdQ3y - cz_half*AdQ3z;


                        // Predict Q at the cell edges (x)
                        // (note that _plus is shifted due to grid location)
                        if ( box_x.contains(i,j,k) ) {
                            Q_minus_x(i,j,k,0) = Q_tilde0 + dQ0x/2.0;
                            Q_minus_x(i,j,k,1) = Q_tilde1 + dQ1x/2.0;
                            Q_minus_x(i,j,k,2) = Q_tilde2 + dQ2x/2.0;
                            Q_minus_x(i,j,k,3) = Q_tilde3 + dQ3x/2.0;
                        }
                        if ( box_x.contains(i-1,j,k) ) {
                            Q_plus_x(i-1,j,k,0) = Q_tilde0 - dQ0x/2.0;
                            Q_plus_x(i-1,j,k,1) = Q_tilde1 - dQ1x/2.0;
                            Q_plus_x(i-1,j,k,2) = Q_tilde2 - dQ2x/2.0;
                            Q_plus_x(i-1,j,k,3) = Q_tilde3 - dQ3x/2.0;
                        }

                        // Positivity and Monotonicty Limiter for density N:
                        if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) )) {
                            if ((Q_minus_x(i,j,k,0) < 0.0) || (Q_plus_x(i-1,j,k,0) < 0.0)){
                                Q_minus_x(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_x(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_x(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_x(i,j,k,3) = NUz_arr(i,j,k);
                                Q_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                                Q_plus_x(i-1,j,k,1) = NUx_arr(i,j,k);
                                Q_plus_x(i-1,j,k,2) = NUy_arr(i,j,k);
                                Q_plus_x(i-1,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) != 1)) {
                            if (Q_minus_x(i,j,k,0) < 0.0) {
                                Q_minus_x(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_x(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_x(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_x(i,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_x.contains(i,j,k) != 1 ) && ( box_x.contains(i-1,j,k) )) {
                            if (Q_plus_x(i-1,j,k,0) < 0.0){
                                Q_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                                Q_plus_x(i-1,j,k,1) = NUx_arr(i,j,k);
                                Q_plus_x(i-1,j,k,2) = NUy_arr(i,j,k);
                                Q_plus_x(i-1,j,k,3) = NUz_arr(i,j,k);
                            }
                        }

                        // Predict Q at the cell edges (y)
                        if ( box_y.contains(i,j,k) ) {
                            Q_minus_y(i,j,k,0) = Q_tilde0 + dQ0y/2.0;
                            Q_minus_y(i,j,k,1) = Q_tilde1 + dQ1y/2.0;
                            Q_minus_y(i,j,k,2) = Q_tilde2 + dQ2y/2.0;
                            Q_minus_y(i,j,k,3) = Q_tilde3 + dQ3y/2.0;
                        }
                        if ( box_y.contains(i,j-1,k) ) {
                            Q_plus_y(i,j-1,k,0) = Q_tilde0 - dQ0y/2.0;
                            Q_plus_y(i,j-1,k,1) = Q_tilde1 - dQ1y/2.0;
                            Q_plus_y(i,j-1,k,2) = Q_tilde2 - dQ2y/2.0;
                            Q_plus_y(i,j-1,k,3) = Q_tilde3 - dQ3y/2.0;
                        }

                        // Positivity and Monotonicty Limiter for density N:
                        if (( box_y.contains(i,j,k) ) && ( box_y.contains(i,j-1,k) )) {
                            if ((Q_minus_y(i,j,k,0) < 0.0) || (Q_plus_y(i,j-1,k,0) < 0.0)){
                                Q_minus_y(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_y(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_y(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_y(i,j,k,3) = NUz_arr(i,j,k);
                                Q_plus_y(i,j-1,k,0) = N_arr(i,j,k);
                                Q_plus_y(i,j-1,k,1) = NUx_arr(i,j,k);
                                Q_plus_y(i,j-1,k,2) = NUy_arr(i,j,k);
                                Q_plus_y(i,j-1,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_y.contains(i,j,k) ) && ( box_y.contains(i,j-1,k) != 1)) {
                            if (Q_minus_y(i,j,k,0) < 0.0) {
                                Q_minus_y(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_y(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_y(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_y(i,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_y.contains(i,j,k) != 1 ) && ( box_y.contains(i,j-1,k) )) {
                            if (Q_plus_y(i,j-1,k,0) < 0.0){
                                Q_plus_y(i,j-1,k,0) = N_arr(i,j,k);
                                Q_plus_y(i,j-1,k,1) = NUx_arr(i,j,k);
                                Q_plus_y(i,j-1,k,2) = NUy_arr(i,j,k);
                                Q_plus_y(i,j-1,k,3) = NUz_arr(i,j,k);
                            }
                        }

                        if ( box_z.contains(i,j,k) ) {
                        // Predict Q at the cell edges (z)
                            Q_minus_z(i,j,k,0) = Q_tilde0 + dQ0z/2.0;
                            Q_minus_z(i,j,k,1) = Q_tilde1 + dQ1z/2.0;
                            Q_minus_z(i,j,k,2) = Q_tilde2 + dQ2z/2.0;
                            Q_minus_z(i,j,k,3) = Q_tilde3 + dQ3z/2.0;
                        }
                        if ( box_z.contains(i,j,k-1) ) {
                            Q_plus_z(i,j,k-1,0) = Q_tilde0 - dQ0z/2.0;
                            Q_plus_z(i,j,k-1,1) = Q_tilde1 - dQ1z/2.0;
                            Q_plus_z(i,j,k-1,2) = Q_tilde2 - dQ2z/2.0;
                            Q_plus_z(i,j,k-1,3) = Q_tilde3 - dQ3z/2.0;
                        }


                        // Positivity and Monotonicty Limiter for density N: z
                        if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j,k-1) )) {
                            if ((Q_minus_z(i,j,k,0) < 0.0) || (Q_plus_z(i,j,k-1,0) < 0.0)){
                                Q_minus_z(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_z(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_z(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_z(i,j,k,3) = NUz_arr(i,j,k);
                                Q_plus_z(i,j,k-1,0) = N_arr(i,j,k);
                                Q_plus_z(i,j,k-1,1) = NUx_arr(i,j,k);
                                Q_plus_z(i,j,k-1,2) = NUy_arr(i,j,k);
                                Q_plus_z(i,j,k-1,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j,k-1) != 1)) {
                            if (Q_minus_z(i,j,k,0) < 0.0) {
                                Q_minus_z(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_z(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_z(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_z(i,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_z.contains(i,j,k) != 1 ) && ( box_z.contains(i,j,k-1) )) {
                            if (Q_plus_z(i,j,k-1,0) < 0.0){
                                Q_plus_z(i,j,k-1,0) = N_arr(i,j,k);
                                Q_plus_z(i,j,k-1,1) = NUx_arr(i,j,k);
                                Q_plus_z(i,j,k-1,2) = NUy_arr(i,j,k);
                                Q_plus_z(i,j,k-1,3) = NUz_arr(i,j,k);
                            }
                        }

                    #elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)

                        // Compute the cell slopes x
                        amrex::Real  dQ0x = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                        amrex::Real  dQ1x = ave( NUx_arr(i,j,k) - NUx_arr(i-1,j,k) , NUx_arr(i+1,j,k) - NUx_arr(i,j,k) );
                        amrex::Real  dQ2x = ave( NUy_arr(i,j,k) - NUy_arr(i-1,j,k) , NUy_arr(i+1,j,k) - NUy_arr(i,j,k) );
                        amrex::Real  dQ3x = ave( NUz_arr(i,j,k) - NUz_arr(i-1,j,k) , NUz_arr(i+1,j,k) - NUz_arr(i,j,k) );

                        // Compute the cell slopes z
                        amrex::Real  dQ0z = ave( N_arr(i,j,k) - N_arr(i,j-1,k) , N_arr(i,j+1,k) - N_arr(i,j,k) );
                        amrex::Real  dQ1z = ave( NUx_arr(i,j,k) - NUx_arr(i,j-1,k) , NUx_arr(i,j+1,k) - NUx_arr(i,j,k) );
                        amrex::Real  dQ2z = ave( NUy_arr(i,j,k) - NUy_arr(i,j-1,k) , NUy_arr(i,j+1,k) - NUy_arr(i,j,k) );
                        amrex::Real  dQ3z = ave( NUz_arr(i,j,k) - NUz_arr(i,j-1,k) , NUz_arr(i,j+1,k) - NUz_arr(i,j,k) );

                        #if defined(WARPX_DIM_RZ)
                        // TODO: Generalize this condition
                        // Impose "none" boundaries
                        // Condition: dQx = 0 at r = 0
                        if  (i == domain.smallEnd(0)) {
                            // TODO BC: Reflected across r = 0:
                            // R|_{0+} -> L|_{0-}
                            // N -> N (N_arr(i-1,j,k) -> N_arr(i+1,j,k))
                            // NUr -> -NUr (NUx_arr(i-1,j,k) -> -NUx_arr(i+1,j,k))
                            // NUt -> -NUt (NUy_arr(i-1,j,k) -> -NUy_arr(i+1,j,k))
                            // NUz -> -NUz (NUz_arr(i-1,j,k) -> NUz_arr(i+1,j,k))
                            dQ0x = ave( N_arr(i,j,k) - N_arr(i+1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                            dQ1x = ave( NUx_arr(i,j,k) + NUx_arr(i+1,j,k) , NUx_arr(i+1,j,k) - NUx_arr(i,j,k) );
                            dQ2x = ave( NUy_arr(i,j,k) + NUy_arr(i+1,j,k) , NUy_arr(i+1,j,k) - NUy_arr(i,j,k) );
                            dQ3x = ave( NUz_arr(i,j,k) - NUz_arr(i+1,j,k) , NUz_arr(i+1,j,k) - NUz_arr(i,j,k) );
                        } else if (i == domain.bigEnd(0)+1) {
                            dQ0x = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , 0.0 );
                            dQ1x = ave( NUx_arr(i,j,k) - NUx_arr(i-1,j,k) , 0.0 );
                            dQ2x = ave( NUy_arr(i,j,k) - NUy_arr(i-1,j,k) , 0.0 );
                            dQ3x = ave( NUz_arr(i,j,k) - NUz_arr(i-1,j,k) , 0.0 );
                        }
                        #endif

                        // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
                        amrex::Real  AdQ0x = A00x*dQ0x + A01x*dQ1x + A02x*dQ2x + A03x*dQ3x;
                        amrex::Real  AdQ1x = A10x*dQ0x + A11x*dQ1x + A12x*dQ2x + A13x*dQ3x;
                        amrex::Real  AdQ2x = A20x*dQ0x + A21x*dQ1x + A22x*dQ2x + A23x*dQ3x;
                        amrex::Real  AdQ3x = A30x*dQ0x + A31x*dQ1x + A32x*dQ2x + A33x*dQ3x;
                        amrex::Real  AdQ0z = A00z*dQ0z + A01z*dQ1z + A02z*dQ2z + A03z*dQ3z;
                        amrex::Real  AdQ1z = A10z*dQ0z + A11z*dQ1z + A12z*dQ2z + A13z*dQ3z;
                        amrex::Real  AdQ2z = A20z*dQ0z + A21z*dQ1z + A22z*dQ2z + A23z*dQ3z;
                        amrex::Real  AdQ3z = A30z*dQ0z + A31z*dQ1z + A32z*dQ2z + A33z*dQ3z;
                        amrex::Real  Q_tilde0 = N_arr(i,j,k)   - cx_half*AdQ0x - cz_half*AdQ0z;
                        amrex::Real  Q_tilde1 = NUx_arr(i,j,k) - cx_half*AdQ1x - cz_half*AdQ1z;
                        amrex::Real  Q_tilde2 = NUy_arr(i,j,k) - cx_half*AdQ2x - cz_half*AdQ2z;
                        amrex::Real  Q_tilde3 = NUz_arr(i,j,k) - cx_half*AdQ3x - cz_half*AdQ3z;

                        // Predict Q at the cell edges (x)
                        // (note that _plus is shifted due to grid location)
                        if ( box_x.contains(i,j,k) ) {
                            Q_minus_x(i,j,k,0) = Q_tilde0 + dQ0x/2.0;
                            Q_minus_x(i,j,k,1) = Q_tilde1 + dQ1x/2.0;
                            Q_minus_x(i,j,k,2) = Q_tilde2 + dQ2x/2.0;
                            Q_minus_x(i,j,k,3) = Q_tilde3 + dQ3x/2.0;
                        }
                        if ( box_x.contains(i-1,j,k) ) {
                            Q_plus_x(i-1,j,k,0) = Q_tilde0 - dQ0x/2.0;
                            Q_plus_x(i-1,j,k,1) = Q_tilde1 - dQ1x/2.0;
                            Q_plus_x(i-1,j,k,2) = Q_tilde2 - dQ2x/2.0;
                            Q_plus_x(i-1,j,k,3) = Q_tilde3 - dQ3x/2.0;
                        }

                    // Positivity and Monotonicty Limiter for density N,
                    // This sets the slope (dQ) to zero for all quantities
                        if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) )) {
                            if ((Q_minus_x(i,j,k,0) < 0.0) || (Q_plus_x(i-1,j,k,0) < 0.0)){
                                Q_minus_x(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_x(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_x(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_x(i,j,k,3) = NUz_arr(i,j,k);
                                Q_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                                Q_plus_x(i-1,j,k,1) = NUx_arr(i,j,k);
                                Q_plus_x(i-1,j,k,2) = NUy_arr(i,j,k);
                                Q_plus_x(i-1,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_x.contains(i,j,k) ) && ( box_x.contains(i-1,j,k) != 1)) {
                            if (Q_minus_x(i,j,k,0) < 0.0) {
                                Q_minus_x(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_x(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_x(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_x(i,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_x.contains(i,j,k) != 1 ) && ( box_x.contains(i-1,j,k) )) {
                            if (Q_plus_x(i-1,j,k,0) < 0.0){
                                Q_plus_x(i-1,j,k,0) = N_arr(i,j,k);
                                Q_plus_x(i-1,j,k,1) = NUx_arr(i,j,k);
                                Q_plus_x(i-1,j,k,2) = NUy_arr(i,j,k);
                                Q_plus_x(i-1,j,k,3) = NUz_arr(i,j,k);
                            }
                        }

                        if ( box_z.contains(i,j,k) ) {
                        // Predict Q at the cell edges (z)
                            Q_minus_z(i,j,k,0) = Q_tilde0 + dQ0z/2.0;
                            Q_minus_z(i,j,k,1) = Q_tilde1 + dQ1z/2.0;
                            Q_minus_z(i,j,k,2) = Q_tilde2 + dQ2z/2.0;
                            Q_minus_z(i,j,k,3) = Q_tilde3 + dQ3z/2.0;
                        }
                        if ( box_z.contains(i,j-1,k) ) {
                            Q_plus_z(i,j-1,k,0) = Q_tilde0 - dQ0z/2.0;
                            Q_plus_z(i,j-1,k,1) = Q_tilde1 - dQ1z/2.0;
                            Q_plus_z(i,j-1,k,2) = Q_tilde2 - dQ2z/2.0;
                            Q_plus_z(i,j-1,k,3) = Q_tilde3 - dQ3z/2.0;
                        }

                        // Positivity and Monotonicty Limiter for density N: z
                        if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j-1,k) )) {
                            if ((Q_minus_z(i,j,k,0) < 0.0) || (Q_plus_z(i,j-1,k,0) < 0.0)){
                                Q_minus_z(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_z(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_z(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_z(i,j,k,3) = NUz_arr(i,j,k);
                                Q_plus_z(i,j-1,k,0) = N_arr(i,j,k);
                                Q_plus_z(i,j-1,k,1) = NUx_arr(i,j,k);
                                Q_plus_z(i,j-1,k,2) = NUy_arr(i,j,k);
                                Q_plus_z(i,j-1,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_z.contains(i,j,k) ) && ( box_z.contains(i,j-1,k) != 1)) {
                            if (Q_minus_z(i,j,k,0) < 0.0) {
                                Q_minus_z(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_z(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_z(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_z(i,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_z.contains(i,j,k) != 1 ) && ( box_z.contains(i,j-1,k) )) {
                            if (Q_plus_z(i,j-1,k,0) < 0.0){
                                Q_plus_z(i,j-1,k,0) = N_arr(i,j,k);
                                Q_plus_z(i,j-1,k,1) = NUx_arr(i,j,k);
                                Q_plus_z(i,j-1,k,2) = NUy_arr(i,j,k);
                                Q_plus_z(i,j-1,k,3) = NUz_arr(i,j,k);
                            }
                        }

                    #else

                        // Compute the cell slopes z
                        amrex::Real  dQ0z = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
                        amrex::Real  dQ1z = ave( NUx_arr(i,j,k) - NUx_arr(i-1,j,k) , NUx_arr(i+1,j,k) - NUx_arr(i,j,k) );
                        amrex::Real  dQ2z = ave( NUy_arr(i,j,k) - NUy_arr(i-1,j,k) , NUy_arr(i+1,j,k) - NUy_arr(i,j,k) );
                        amrex::Real  dQ3z = ave( NUz_arr(i,j,k) - NUz_arr(i-1,j,k) , NUz_arr(i+1,j,k) - NUz_arr(i,j,k) );

                        // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
                        amrex::Real  AdQ0z = A00z*dQ0z + A01z*dQ1z + A02z*dQ2z + A03z*dQ3z;
                        amrex::Real  AdQ1z = A10z*dQ0z + A11z*dQ1z + A12z*dQ2z + A13z*dQ3z;
                        amrex::Real  AdQ2z = A20z*dQ0z + A21z*dQ1z + A22z*dQ2z + A23z*dQ3z;
                        amrex::Real  AdQ3z = A30z*dQ0z + A31z*dQ1z + A32z*dQ2z + A33z*dQ3z;
                        amrex::Real  Q_tilde0 = N_arr(i,j,k)   - cz_half*AdQ0z;
                        amrex::Real  Q_tilde1 = NUx_arr(i,j,k) - cz_half*AdQ1z;
                        amrex::Real  Q_tilde2 = NUy_arr(i,j,k) - cz_half*AdQ2z;
                        amrex::Real  Q_tilde3 = NUz_arr(i,j,k) - cz_half*AdQ3z;

                        // Predict Q at the cell edges (z)
                        // (note that _plus is shifted due to grid location)
                        if ( box_z.contains(i,j,k) ) {
                        // Predict Q at the cell edges (z)
                            Q_minus_z(i,j,k,0) = Q_tilde0 + dQ0z/2.0;
                            Q_minus_z(i,j,k,1) = Q_tilde1 + dQ1z/2.0;
                            Q_minus_z(i,j,k,2) = Q_tilde2 + dQ2z/2.0;
                            Q_minus_z(i,j,k,3) = Q_tilde3 + dQ3z/2.0;
                        }
                        if ( box_z.contains(i-1,j,k) ) {
                            Q_plus_z(i-1,j,k,0) = Q_tilde0 - dQ0z/2.0;
                            Q_plus_z(i-1,j,k,1) = Q_tilde1 - dQ1z/2.0;
                            Q_plus_z(i-1,j,k,2) = Q_tilde2 - dQ2z/2.0;
                            Q_plus_z(i-1,j,k,3) = Q_tilde3 - dQ3z/2.0;
                        }

                        // Positivity and Monotonicty Limiter for density N,
                        // This sets the slope (dQ) to zero for all quantities
                        if (( box_z.contains(i,j,k) ) && ( box_z.contains(i-1,j,k) )) {
                            if ((Q_minus_z(i,j,k,0) < 0.0) || (Q_plus_z(i-1,j,k,0) < 0.0)) {
                                Q_minus_z(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_z(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_z(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_z(i,j,k,3) = NUz_arr(i,j,k);
                                Q_plus_z(i-1,j,k,0) = N_arr(i,j,k);
                                Q_plus_z(i-1,j,k,1) = NUx_arr(i,j,k);
                                Q_plus_z(i-1,j,k,2) = NUy_arr(i,j,k);
                                Q_plus_z(i-1,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_z.contains(i,j,k) ) && ( box_z.contains(i-1,j,k) != 1)) {
                            if (Q_minus_z(i,j,k,0) < 0.0) {
                                Q_minus_z(i,j,k,0) = N_arr(i,j,k);
                                Q_minus_z(i,j,k,1) = NUx_arr(i,j,k);
                                Q_minus_z(i,j,k,2) = NUy_arr(i,j,k);
                                Q_minus_z(i,j,k,3) = NUz_arr(i,j,k);
                            }
                        } else if (( box_z.contains(i,j,k) != 1 ) && ( box_z.contains(i-1,j,k) )) {
                            if (Q_plus_z(i-1,j,k,0) < 0.0){
                                Q_plus_z(i-1,j,k,0) = N_arr(i,j,k);
                                Q_plus_z(i-1,j,k,1) = NUx_arr(i,j,k);
                                Q_plus_z(i-1,j,k,2) = NUy_arr(i,j,k);
                                Q_plus_z(i-1,j,k,3) = NUz_arr(i,j,k);
                            }
                        }

                    #endif
                // If N<= 0 then set the boundaries to zero
                } else {
                    #if defined(WARPX_DIM_3D) // 3D:
                    if ( box_x.contains(i,j,k) ) {
                        Q_minus_x(i,j,k,0) = 0.0;
                        Q_minus_x(i,j,k,1) = 0.0;
                        Q_minus_x(i,j,k,2) = 0.0;
                        Q_minus_x(i,j,k,3) = 0.0;
                    }
                    if ( box_x.contains(i-1,j,k) ) {
                        Q_plus_x(i-1,j,k,0) = 0.0;
                        Q_plus_x(i-1,j,k,1) = 0.0;
                        Q_plus_x(i-1,j,k,2) = 0.0;
                        Q_plus_x(i-1,j,k,3) = 0.0;
                    }
                    if ( box_y.contains(i,j,k) ) {
                        Q_minus_y(i,j,k,0) = 0.0;
                        Q_minus_y(i,j,k,1) = 0.0;
                        Q_minus_y(i,j,k,2) = 0.0;
                        Q_minus_y(i,j,k,3) = 0.0;
                    }
                    if ( box_y.contains(i,j-1,k) ) {
                        Q_plus_y(i,j-1,k,0) = 0.0;
                        Q_plus_y(i,j-1,k,1) = 0.0;
                        Q_plus_y(i,j-1,k,2) = 0.0;
                        Q_plus_y(i,j-1,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j,k) ) {
                        Q_minus_z(i,j,k,0) = 0.0;
                        Q_minus_z(i,j,k,1) = 0.0;
                        Q_minus_z(i,j,k,2) = 0.0;
                        Q_minus_z(i,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j,k-1) ) {
                        Q_plus_z(i,j,k-1,0) = 0.0;
                        Q_plus_z(i,j,k-1,1) = 0.0;
                        Q_plus_z(i,j,k-1,2) = 0.0;
                        Q_plus_z(i,j,k-1,3) = 0.0;
                    }
                    #elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ) // 2D:
                    if ( box_x.contains(i,j,k) ) {
                        Q_minus_x(i,j,k,0) = 0.0;
                        Q_minus_x(i,j,k,1) = 0.0;
                        Q_minus_x(i,j,k,2) = 0.0;
                        Q_minus_x(i,j,k,3) = 0.0;
                    }
                    if ( box_x.contains(i-1,j,k) ) {
                        Q_plus_x(i-1,j,k,0) = 0.0;
                        Q_plus_x(i-1,j,k,1) = 0.0;
                        Q_plus_x(i-1,j,k,2) = 0.0;
                        Q_plus_x(i-1,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j,k) ) {
                        Q_minus_z(i,j,k,0) = 0.0;
                        Q_minus_z(i,j,k,1) = 0.0;
                        Q_minus_z(i,j,k,2) = 0.0;
                        Q_minus_z(i,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i,j-1,k) ) {
                        Q_plus_z(i,j-1,k,0) = 0.0;
                        Q_plus_z(i,j-1,k,1) = 0.0;
                        Q_plus_z(i,j-1,k,2) = 0.0;
                        Q_plus_z(i,j-1,k,3) = 0.0;
                    }
                    #else // 1D:
                    if ( box_z.contains(i,j,k) ) {
                        Q_minus_z(i,j,k,0) = 0.0;
                        Q_minus_z(i,j,k,1) = 0.0;
                        Q_minus_z(i,j,k,2) = 0.0;
                        Q_minus_z(i,j,k,3) = 0.0;
                    }
                    if ( box_z.contains(i-1,j,k) ) {
                        Q_plus_z(i-1,j,k,0) = 0.0;
                        Q_plus_z(i-1,j,k,1) = 0.0;
                        Q_plus_z(i-1,j,k,2) = 0.0;
                        Q_plus_z(i-1,j,k,3) = 0.0;
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
            amrex::Array4<amrex::Real> const &Q_minus_x = tmp_Q_minus_x.array(mfi);
            amrex::Array4<amrex::Real> const &Q_plus_x = tmp_Q_plus_x.array(mfi);
            amrex::Array4<amrex::Real> const &Q_minus_y = tmp_Q_minus_y.array(mfi);
            amrex::Array4<amrex::Real> const &Q_plus_y = tmp_Q_plus_y.array(mfi);
            amrex::Array4<amrex::Real> const &Q_minus_z = tmp_Q_minus_z.array(mfi);
            amrex::Array4<amrex::Real> const &Q_plus_z = tmp_Q_plus_z.array(mfi);
        #elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            amrex::Array4<amrex::Real> const &Q_minus_x = tmp_Q_minus_x.array(mfi);
            amrex::Array4<amrex::Real> const &Q_plus_x = tmp_Q_plus_x.array(mfi);
            amrex::Array4<amrex::Real> const &Q_minus_z = tmp_Q_minus_z.array(mfi);
            amrex::Array4<amrex::Real> const &Q_plus_z = tmp_Q_plus_z.array(mfi);
        #else
            amrex::Array4<amrex::Real> const &Q_minus_z = tmp_Q_minus_z.array(mfi);
            amrex::Array4<amrex::Real> const &Q_plus_z = tmp_Q_plus_z.array(mfi);
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
                    if (Q_minus_x(i-1,j,k,0)>0.0) Vx_L_minus = V_calc(Q_minus_x(i-1,j,k,0),Q_minus_x(i-1,j,k,1),Q_minus_x(i-1,j,k,2),Q_minus_x(i-1,j,k,3),clight,0);
                    if (Q_minus_x(i,j,k,0)>0.0)   Vx_I_minus = V_calc(Q_minus_x(i,j,k,0),Q_minus_x(i,j,k,1),Q_minus_x(i,j,k,2),Q_minus_x(i,j,k,3),clight,0);
                    if (Q_plus_x(i-1,j,k,0)>0.0)   Vx_L_plus = V_calc(Q_plus_x(i-1,j,k,0),Q_plus_x(i-1,j,k,1),Q_plus_x(i-1,j,k,2),Q_plus_x(i-1,j,k,3),clight,0);
                    if (Q_plus_x(i,j,k,0)>0.0) Vx_I_plus = V_calc(Q_plus_x(i,j,k,0),Q_plus_x(i,j,k,1),Q_plus_x(i,j,k,2),Q_plus_x(i,j,k,3),clight,0);

                    if (Q_minus_y(i,j-1,k,0)>0.0) Vy_L_minus = V_calc(Q_minus_y(i,j-1,k,0),Q_minus_y(i,j-1,k,1),Q_minus_y(i,j-1,k,2),Q_minus_y(i,j-1,k,3),clight,1);
                    if (Q_minus_y(i,j,k,0)>0.0)   Vy_I_minus = V_calc(Q_minus_y(i,j,k,0),Q_minus_y(i,j,k,1),Q_minus_y(i,j,k,2),Q_minus_y(i,j,k,3),clight,1);
                    if (Q_plus_y(i,j-1,k,0)>0.0)   Vy_L_plus = V_calc(Q_plus_y(i,j-1,k,0),Q_plus_y(i,j-1,k,1),Q_plus_y(i,j-1,k,2),Q_plus_y(i,j-1,k,3),clight,1);
                    if (Q_plus_y(i,j,k,0)>0.0)Vy_I_plus = V_calc(Q_plus_y(i,j,k,0),Q_plus_y(i,j,k,1),Q_plus_y(i,j,k,2),Q_plus_y(i,j,k,3),clight,1);

                    if (Q_minus_z(i,j,k-1,0)>0.0) Vz_L_minus = V_calc(Q_minus_z(i,j,k-1,0),Q_minus_z(i,j,k-1,1),Q_minus_z(i,j,k-1,2),Q_minus_z(i,j,k-1,3),clight,2);
                    if (Q_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(Q_minus_z(i,j,k,0),Q_minus_z(i,j,k,1),Q_minus_z(i,j,k,2),Q_minus_z(i,j,k,3),clight,2);
                    if (Q_plus_z(i,j,k-1,0)>0.0)   Vz_L_plus = V_calc(Q_plus_z(i,j,k-1,0),Q_plus_z(i,j,k-1,1),Q_plus_z(i,j,k-1,2),Q_plus_z(i,j,k-1,3),clight,2);
                    if (Q_plus_z(i,j,k,0)>0.0) Vz_I_plus = V_calc(Q_plus_z(i,j,k,0),Q_plus_z(i,j,k,1),Q_plus_z(i,j,k,2),Q_plus_z(i,j,k,3),clight,2);

                    // compute the fluxes:
                    // (note that _plus is shifted due to grid location)
                    amrex::Real F0_minusx = flux(Q_minus_x(i-1,j,k,0),Q_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F0_plusx =  flux(Q_minus_x(i,j,k,0),  Q_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                    amrex::Real F1_minusx = flux(Q_minus_x(i-1,j,k,1),Q_plus_x(i-1,j,k,1),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F1_plusx =  flux(Q_minus_x(i,j,k,1),  Q_plus_x(i,j,k,1),    Vx_I_minus,Vx_I_plus);
                    amrex::Real F2_minusx = flux(Q_minus_x(i-1,j,k,2),Q_plus_x(i-1,j,k,2),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F2_plusx =  flux(Q_minus_x(i,j,k,2),  Q_plus_x(i,j,k,2),    Vx_I_minus,Vx_I_plus);
                    amrex::Real F3_minusx = flux(Q_minus_x(i-1,j,k,3),Q_plus_x(i-1,j,k,3),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F3_plusx =  flux(Q_minus_x(i,j,k,3),  Q_plus_x(i,j,k,3),    Vx_I_minus,Vx_I_plus);

                    amrex::Real F0_minusy = flux(Q_minus_y(i,j-1,k,0),Q_plus_y(i,j-1,k,0),  Vy_L_minus,Vy_L_plus);
                    amrex::Real F0_plusy =  flux(Q_minus_y(i,j,k,0),  Q_plus_y(i,j,k,0),    Vy_I_minus,Vy_I_plus);
                    amrex::Real F1_minusy = flux(Q_minus_y(i,j-1,k,1),Q_plus_y(i,j-1,k,1),  Vy_L_minus,Vy_L_plus);
                    amrex::Real F1_plusy =  flux(Q_minus_y(i,j,k,1),  Q_plus_y(i,j,k,1),    Vy_I_minus,Vy_I_plus);
                    amrex::Real F2_minusy = flux(Q_minus_y(i,j-1,k,2),Q_plus_y(i,j-1,k,2),  Vy_L_minus,Vy_L_plus);
                    amrex::Real F2_plusy =  flux(Q_minus_y(i,j,k,2),  Q_plus_y(i,j,k,2),    Vy_I_minus,Vy_I_plus);
                    amrex::Real F3_minusy = flux(Q_minus_y(i,j-1,k,3),Q_plus_y(i,j-1,k,3),  Vy_L_minus,Vy_L_plus);
                    amrex::Real F3_plusy =  flux(Q_minus_y(i,j,k,3),  Q_plus_y(i,j,k,3),    Vy_I_minus,Vy_I_plus);

                    amrex::Real F0_minusz = flux(Q_minus_z(i,j,k-1,0),Q_plus_z(i,j,k-1,0),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F0_plusz =  flux(Q_minus_z(i,j,k,0),  Q_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F1_minusz = flux(Q_minus_z(i,j,k-1,1),Q_plus_z(i,j,k-1,1),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F1_plusz =  flux(Q_minus_z(i,j,k,1),  Q_plus_z(i,j,k,1),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F2_minusz = flux(Q_minus_z(i,j,k-1,2),Q_plus_z(i,j,k-1,2),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F2_plusz =  flux(Q_minus_z(i,j,k,2),  Q_plus_z(i,j,k,2),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F3_minusz = flux(Q_minus_z(i,j,k-1,3),Q_plus_z(i,j,k-1,3),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F3_plusz =  flux(Q_minus_z(i,j,k,3),  Q_plus_z(i,j,k,3),    Vz_I_minus,Vz_I_plus);

                    // Update Q from tn -> tn + dt
                    N_arr(i,j,k) = N_arr(i,j,k) - cx*(F0_plusx - F0_minusx)
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

                #elif defined(WARPX_DIM_XZ)

                    amrex::Real Vx_L_minus = 0.0, Vx_I_minus = 0.0, Vx_L_plus = 0.0, Vx_I_plus = 0.0;
                    amrex::Real Vz_L_minus = 0.0, Vz_I_minus = 0.0, Vz_L_plus = 0.0, Vz_I_plus = 0.0;

                    // Verify positive density, then compute velocity
                    if (Q_minus_x(i-1,j,k,0)>0.0) Vx_L_minus = V_calc(Q_minus_x(i-1,j,k,0),Q_minus_x(i-1,j,k,1),Q_minus_x(i-1,j,k,2),Q_minus_x(i-1,j,k,3),clight,0);
                    if (Q_minus_x(i,j,k,0)>0.0)   Vx_I_minus = V_calc(Q_minus_x(i,j,k,0),Q_minus_x(i,j,k,1),Q_minus_x(i,j,k,2),Q_minus_x(i,j,k,3),clight,0);
                    if (Q_plus_x(i-1,j,k,0)>0.0)   Vx_L_plus = V_calc(Q_plus_x(i-1,j,k,0),Q_plus_x(i-1,j,k,1),Q_plus_x(i-1,j,k,2),Q_plus_x(i-1,j,k,3),clight,0);
                    if (Q_plus_x(i,j,k,0)>0.0) Vx_I_plus = V_calc(Q_plus_x(i,j,k,0),Q_plus_x(i,j,k,1),Q_plus_x(i,j,k,2),Q_plus_x(i,j,k,3),clight,0);

                    if (Q_minus_z(i,j-1,k,0)>0.0) Vz_L_minus = V_calc(Q_minus_z(i,j-1,k,0),Q_minus_z(i,j-1,k,1),Q_minus_z(i,j-1,k,2),Q_minus_z(i,j-1,k,3),clight,2);
                    if (Q_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(Q_minus_z(i,j,k,0),Q_minus_z(i,j,k,1),Q_minus_z(i,j,k,2),Q_minus_z(i,j,k,3),clight,2);
                    if (Q_plus_z(i,j-1,k,0)>0.0)   Vz_L_plus = V_calc(Q_plus_z(i,j-1,k,0),Q_plus_z(i,j-1,k,1),Q_plus_z(i,j-1,k,2),Q_plus_z(i,j-1,k,3),clight,2);
                    if (Q_plus_z(i,j,k,0)>0.0) Vz_I_plus = V_calc(Q_plus_z(i,j,k,0),Q_plus_z(i,j,k,1),Q_plus_z(i,j,k,2),Q_plus_z(i,j,k,3),clight,2);


                    // compute the fluxes:
                    // (note that _plus is shifted due to grid location)
                    amrex::Real F0_minusx = flux(Q_minus_x(i-1,j,k,0),Q_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F0_plusx =  flux(Q_minus_x(i,j,k,0),  Q_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus);
                    amrex::Real F1_minusx = flux(Q_minus_x(i-1,j,k,1),Q_plus_x(i-1,j,k,1),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F1_plusx =  flux(Q_minus_x(i,j,k,1),  Q_plus_x(i,j,k,1),    Vx_I_minus,Vx_I_plus);
                    amrex::Real F2_minusx = flux(Q_minus_x(i-1,j,k,2),Q_plus_x(i-1,j,k,2),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F2_plusx =  flux(Q_minus_x(i,j,k,2),  Q_plus_x(i,j,k,2),    Vx_I_minus,Vx_I_plus);
                    amrex::Real F3_minusx = flux(Q_minus_x(i-1,j,k,3),Q_plus_x(i-1,j,k,3),  Vx_L_minus,Vx_L_plus);
                    amrex::Real F3_plusx =  flux(Q_minus_x(i,j,k,3),  Q_plus_x(i,j,k,3),    Vx_I_minus,Vx_I_plus);

                    amrex::Real F0_minusz = flux(Q_minus_z(i,j-1,k,0),Q_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F0_plusz =  flux(Q_minus_z(i,j,k,0),  Q_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F1_minusz = flux(Q_minus_z(i,j-1,k,1),Q_plus_z(i,j-1,k,1),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F1_plusz =  flux(Q_minus_z(i,j,k,1),  Q_plus_z(i,j,k,1),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F2_minusz = flux(Q_minus_z(i,j-1,k,2),Q_plus_z(i,j-1,k,2),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F2_plusz =  flux(Q_minus_z(i,j,k,2),  Q_plus_z(i,j,k,2),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F3_minusz = flux(Q_minus_z(i,j-1,k,3),Q_plus_z(i,j-1,k,3),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F3_plusz =  flux(Q_minus_z(i,j,k,3),  Q_plus_z(i,j,k,3),    Vz_I_minus,Vz_I_plus);

                    // Update Q from tn -> tn + dt
                    N_arr(i,j,k) = N_arr(i,j,k) - cx*(F0_plusx - F0_minusx)
                                                - cz*(F0_plusz - F0_minusz);
                    NUx_arr(i,j,k) = NUx_arr(i,j,k) - cx*(F1_plusx - F1_minusx)
                                                    - cz*(F1_plusz - F1_minusz);
                    NUy_arr(i,j,k) = NUy_arr(i,j,k) - cx*(F2_plusx - F2_minusx)
                                                    - cz*(F2_plusz - F2_minusz);
                    NUz_arr(i,j,k) = NUz_arr(i,j,k) - cx*(F3_plusx - F3_minusx)
                                                    - cz*(F3_plusz - F3_minusz);

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
                    } else if (i == domain.bigEnd(0)+1) { // TODO: Fix domain? off by one
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

                    // TODO: Generalize this condition
                    // Impose "none" boundaries
                    // Condition: Vx(r) = 0 at boundaries
                    amrex::Real Vx_L_minus = 0.0, Vx_I_minus = 0.0, Vx_L_plus = 0.0, Vx_I_plus = 0.0;
                    amrex::Real Vz_L_minus = 0.0, Vz_I_minus = 0.0, Vz_L_plus = 0.0, Vz_I_plus = 0.0;
                    if (Q_minus_x(i,j,k,0)>0.0) Vx_I_minus = V_calc(Q_minus_x(i,j,k,0),Q_minus_x(i,j,k,1),Q_minus_x(i,j,k,2),Q_minus_x(i,j,k,3),clight,0);
                    if (Q_plus_x(i-1,j,k,0)>0.0) Vx_L_plus = V_calc(Q_plus_x(i-1,j,k,0),Q_plus_x(i-1,j,k,1),Q_plus_x(i-1,j,k,2),Q_plus_x(i-1,j,k,3),clight,0);

                    if (Q_minus_z(i,j-1,k,0)>0.0) Vz_L_minus = V_calc(Q_minus_z(i,j-1,k,0),Q_minus_z(i,j-1,k,1),Q_minus_z(i,j-1,k,2),Q_minus_z(i,j-1,k,3),clight,2);
                    if (Q_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(Q_minus_z(i,j,k,0),Q_minus_z(i,j,k,1),Q_minus_z(i,j,k,2),Q_minus_z(i,j,k,3),clight,2);
                    if (Q_plus_z(i,j-1,k,0)>0.0)   Vz_L_plus = V_calc(Q_plus_z(i,j-1,k,0),Q_plus_z(i,j-1,k,1),Q_plus_z(i,j-1,k,2),Q_plus_z(i,j-1,k,3),clight,2);
                    if (Q_plus_z(i,j,k,0)>0.0) Vz_I_plus = V_calc(Q_plus_z(i,j,k,0),Q_plus_z(i,j,k,1),Q_plus_z(i,j,k,2),Q_plus_z(i,j,k,3),clight,2);


                    // compute the fluxes:
                    // (note that _plus is shifted due to grid location)
                    amrex::Real F0_minusx = 0.0, F1_minusx = 0.0, F2_minusx = 0.0, F3_minusx = 0.0;
                    amrex::Real F0_plusx = 0.0, F1_plusx = 0.0, F2_plusx = 0.0, F3_plusx = 0.0;
                    if (i != domain.smallEnd(0)) {
                        if (Q_minus_x(i-1,j,k,0)>0.0) Vx_L_minus = V_calc(Q_minus_x(i-1,j,k,0),Q_minus_x(i-1,j,k,1),Q_minus_x(i-1,j,k,2),Q_minus_x(i-1,j,k,3),clight,0);
                        F0_minusx = flux(Q_minus_x(i-1,j,k,0),Q_plus_x(i-1,j,k,0),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                        F1_minusx = flux(Q_minus_x(i-1,j,k,1),Q_plus_x(i-1,j,k,1),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                        F2_minusx = flux(Q_minus_x(i-1,j,k,2),Q_plus_x(i-1,j,k,2),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                        F3_minusx = flux(Q_minus_x(i-1,j,k,3),Q_plus_x(i-1,j,k,3),  Vx_L_minus,Vx_L_plus)*S_Ar_minus;
                    }
                    if (i < domain.bigEnd(0)) {
                        if (Q_plus_x(i,j,k,0)>0.0) Vx_I_plus = V_calc(Q_plus_x(i,j,k,0),Q_plus_x(i,j,k,1),Q_plus_x(i,j,k,2),Q_plus_x(i,j,k,3),clight,0);
                        F0_plusx =  flux(Q_minus_x(i,j,k,0),  Q_plus_x(i,j,k,0),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                        F1_plusx =  flux(Q_minus_x(i,j,k,1),  Q_plus_x(i,j,k,1),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                        F2_plusx =  flux(Q_minus_x(i,j,k,2),  Q_plus_x(i,j,k,2),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                        F3_plusx =  flux(Q_minus_x(i,j,k,3),  Q_plus_x(i,j,k,3),    Vx_I_minus,Vx_I_plus)*S_Ar_plus;
                    }

                    amrex::Real F0_minusz = flux(Q_minus_z(i,j-1,k,0),Q_plus_z(i,j-1,k,0),  Vz_L_minus,Vz_L_plus)*S_Az;
                    amrex::Real F0_plusz =  flux(Q_minus_z(i,j,k,0),  Q_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus)*S_Az;
                    amrex::Real F1_minusz = flux(Q_minus_z(i,j-1,k,1),Q_plus_z(i,j-1,k,1),  Vz_L_minus,Vz_L_plus)*S_Az;
                    amrex::Real F1_plusz =  flux(Q_minus_z(i,j,k,1),  Q_plus_z(i,j,k,1),    Vz_I_minus,Vz_I_plus)*S_Az;
                    amrex::Real F2_minusz = flux(Q_minus_z(i,j-1,k,2),Q_plus_z(i,j-1,k,2),  Vz_L_minus,Vz_L_plus)*S_Az;
                    amrex::Real F2_plusz =  flux(Q_minus_z(i,j,k,2),  Q_plus_z(i,j,k,2),    Vz_I_minus,Vz_I_plus)*S_Az;
                    amrex::Real F3_minusz = flux(Q_minus_z(i,j-1,k,3),Q_plus_z(i,j-1,k,3),  Vz_L_minus,Vz_L_plus)*S_Az;
                    amrex::Real F3_plusz =  flux(Q_minus_z(i,j,k,3),  Q_plus_z(i,j,k,3),    Vz_I_minus,Vz_I_plus)*S_Az;

                    // Update Q from tn -> tn + dt
                    N_arr(i,j,k) = N_arr(i,j,k)     - (dt/Vij)*(F0_plusx - F0_minusx + F0_plusz - F0_minusz);
                    NUx_arr(i,j,k) = NUx_arr(i,j,k) - (dt/Vij)*(F1_plusx - F1_minusx + F1_plusz - F1_minusz);
                    NUy_arr(i,j,k) = NUy_arr(i,j,k) - (dt/Vij)*(F2_plusx - F2_minusx + F2_plusz - F2_minusz);
                    NUz_arr(i,j,k) = NUz_arr(i,j,k) - (dt/Vij)*(F3_plusx - F3_minusx + F3_plusz - F3_minusz);

                #else

                    amrex::Real Vz_L_minus = 0.0, Vz_I_minus = 0.0, Vz_L_plus = 0.0, Vz_I_plus = 0.0;

                    // Compute the half-timestep velocities
                    if (Q_minus_z(i-1,j,k,0)>0.0) Vz_L_minus = V_calc(Q_minus_z(i-1,j,k,0),Q_minus_z(i-1,j,k,1),Q_minus_z(i-1,j,k,2),Q_minus_z(i-1,j,k,3),clight,2);
                    if (Q_minus_z(i,j,k,0)>0.0)   Vz_I_minus = V_calc(Q_minus_z(i,j,k,0),Q_minus_z(i,j,k,1),Q_minus_z(i,j,k,2),Q_minus_z(i,j,k,3),clight,2);
                    if (Q_plus_z(i-1,j,k,0)>0.0)   Vz_L_plus = V_calc(Q_plus_z(i-1,j,k,0),Q_plus_z(i-1,j,k,1),Q_plus_z(i-1,j,k,2),Q_plus_z(i-1,j,k,3),clight,2);
                    if (Q_plus_z(i,j,k,0)>0.0) Vz_I_plus = V_calc(Q_plus_z(i,j,k,0),Q_plus_z(i,j,k,1),Q_plus_z(i,j,k,2),Q_plus_z(i,j,k,3),clight,2);

                    // compute the fluzes:
                    // (note that _plus is shifted due to grid location)
                    amrex::Real F0_minusz = flux(Q_minus_z(i-1,j,k,0),Q_plus_z(i-1,j,k,0),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F0_plusz =  flux(Q_minus_z(i,j,k,0),  Q_plus_z(i,j,k,0),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F1_minusz = flux(Q_minus_z(i-1,j,k,1),Q_plus_z(i-1,j,k,1),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F1_plusz =  flux(Q_minus_z(i,j,k,1),  Q_plus_z(i,j,k,1),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F2_minusz = flux(Q_minus_z(i-1,j,k,2),Q_plus_z(i-1,j,k,2),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F2_plusz =  flux(Q_minus_z(i,j,k,2),  Q_plus_z(i,j,k,2),    Vz_I_minus,Vz_I_plus);
                    amrex::Real F3_minusz = flux(Q_minus_z(i-1,j,k,3),Q_plus_z(i-1,j,k,3),  Vz_L_minus,Vz_L_plus);
                    amrex::Real F3_plusz =  flux(Q_minus_z(i,j,k,3),  Q_plus_z(i,j,k,3),    Vz_I_minus,Vz_I_plus);


                    // Update Q from tn -> tn + dt
                    N_arr(i,j,k) = N_arr(i,j,k)     - cz*(F0_plusz - F0_minusz);
                    NUx_arr(i,j,k) = NUx_arr(i,j,k) - cz*(F1_plusz - F1_minusz);
                    NUy_arr(i,j,k) = NUy_arr(i,j,k) - cz*(F2_plusz - F2_minusz);
                    NUz_arr(i,j,k) = NUz_arr(i,j,k) - cz*(F3_plusz - F3_minusz);
                #endif
            }
        );
    }
}


// Momentum source due to curvature
void WarpXFluidContainer::centrifugal_source (int lev)
{
    WARPX_PROFILE("WarpXFluidContainer::centrifugal_source");

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
                    auto u_r =     (NUx_arr(i, j, k) / (N_arr(i,j,k) * clight ));
                    auto u_theta = (NUy_arr(i, j, k) / (N_arr(i,j,k) * clight ));
                    auto u_z =     (NUz_arr(i, j, k) / (N_arr(i,j,k) * clight ));

                    // (SSP-RK3) Push the fluid momentum (R and Theta)
                    // F_r, F_theta are first order euler pushes of our rhs operator
                    // TODO: only do this if (r != 0)
                    if (i != domain.smallEnd(0)) {
                        auto u_r_1     = F_r(r,u_r,u_theta,u_z,dt);
                        auto u_theta_1 = F_theta(r,u_r,u_theta,u_z,dt);
                        auto u_r_2     = (0.75)*(u_r)     + (0.25)*F_r(r,u_r_1,u_theta_1,u_z,dt);
                        auto u_theta_2 = (0.75)*(u_theta) + (0.25)*F_theta(r,u_r_1,u_theta_1,u_z,dt);
                        u_r            = (1.0/3.0)*(u_r)     + (2.0/3.0)*F_r(r,u_r_2,u_theta_2,u_z,dt);
                        u_theta        = (1.0/3.0)*(u_theta) + (2.0/3.0)*F_theta(r,u_r_2,u_theta_2,u_z,dt);

                        // Calculate NU, save NUr, NUtheta
                        NUx_arr(i,j,k) = N_arr(i,j,k)*u_r*clight;
                        NUy_arr(i,j,k) = N_arr(i,j,k)*u_theta*clight;

                    // TODO FIX: BC r = 0, u_theta = 0, and there is no extra source terms
                    } else {
                        NUx_arr(i,j,k) = 0.0;
                        NUy_arr(i,j,k) = 0.0;
                    }
                }
            }
        );
    }
}

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
    //Check whether m_E_ext_s is "none"
    bool external_e_fields; // Needs intializing
    bool external_b_fields; // Needs intializing


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

    // External field parsers
    external_e_fields = (m_E_ext_s == "parse_e_ext_function");
    external_b_fields = (m_B_ext_s == "parse_b_ext_function");
    if (external_e_fields){
        constexpr auto num_arguments = 4; //x,y,z,t
        m_Exfield_parser = m_Ex_parser->compile<num_arguments>();
        m_Eyfield_parser = m_Ey_parser->compile<num_arguments>();
        m_Ezfield_parser = m_Ez_parser->compile<num_arguments>();
    }

    if (external_b_fields){
        constexpr auto num_arguments = 4; //x,y,z,t
        m_Bxfield_parser = m_Bx_parser->compile<num_arguments>();
        m_Byfield_parser = m_By_parser->compile<num_arguments>();
        m_Bzfield_parser = m_Bz_parser->compile<num_arguments>();
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

                    if (WarpX::gamma_boost > 1._rt) { // Lorentz transform fields due to moving frame
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
                            amrex::Real t_lab = WarpX::gamma_boost*(t + WarpX::beta_boost*z/PhysConst::c);
                            amrex::Real z_lab = WarpX::gamma_boost*(z + WarpX::beta_boost*PhysConst::c*t);

                            // Grab the external fields in the lab frame:
                            if ( external_e_fields ) {
                                Ex_ext_lab = m_Exfield_parser(x, y, z_lab, t_lab);
                                Ey_ext_lab = m_Eyfield_parser(x, y, z_lab, t_lab);
                                Ez_ext_lab = m_Ezfield_parser(x, y, z_lab, t_lab);
                            }else{
                                Ex_ext_lab = 0.0;
                                Ey_ext_lab = 0.0;
                                Ez_ext_lab = 0.0;
                            }
                            if ( external_b_fields ) {
                                Bx_ext_lab = m_Bxfield_parser(x, y, z_lab, t_lab);
                                By_ext_lab = m_Byfield_parser(x, y, z_lab, t_lab);
                                Bz_ext_lab = m_Bzfield_parser(x, y, z_lab, t_lab);
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
                            Ex_ext_boost = WarpX::gamma_boost*(Ex_ext_lab - WarpX::beta_boost*PhysConst::c*By_ext_lab);
                            Ey_ext_boost = WarpX::gamma_boost*(Ey_ext_lab + WarpX::beta_boost*PhysConst::c*Bx_ext_lab);
                            Bx_ext_boost = WarpX::gamma_boost*(Bx_ext_lab + WarpX::beta_boost*Ey_ext_lab/PhysConst::c);
                            By_ext_boost = WarpX::gamma_boost*(By_ext_lab - WarpX::beta_boost*Ex_ext_lab/PhysConst::c);

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

                            Ex_Nodal += m_Exfield_parser(x, y, z, t);
                            Ey_Nodal += m_Eyfield_parser(x, y, z, t);
                            Ez_Nodal += m_Ezfield_parser(x, y, z, t);
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

                            Bx_Nodal += m_Bxfield_parser(x, y, z, t);
                            By_Nodal += m_Byfield_parser(x, y, z, t);
                            Bz_Nodal += m_Bzfield_parser(x, y, z, t);
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
        auto owner_mask_rho_arr = owner_mask_rho->array(mfi);

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
                amrex::Real gamma = 1.0;
                if (N_arr(i, j, k)>0.0) gamma = std::sqrt(N_arr(i, j, k) * N_arr(i, j, k) + (NUx_arr(i, j, k) * NUx_arr(i, j, k) + NUy_arr(i, j, k) * NUy_arr(i, j, k) + NUz_arr(i, j, k) * NUz_arr(i, j, k)) * inv_clight_sq) / N_arr(i, j, k);
                // If density is too small, the result can be undefined, so we need to correct gamma
                if (gamma < 1.0) gamma = 1.0;
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
