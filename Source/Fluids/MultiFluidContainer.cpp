/* Copyright 2023 Grant Johnson, Remi Lehe
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Utils/Parser/ParserUtils.H"

#include <string>

using namespace amrex;

MultiFluidContainer::MultiFluidContainer (int nlevs_max, const amrex::Geometry& geom)
{
    ReadParameters();

    const int nspecies = static_cast<int>(species_names.size());

    allcontainers.resize(nspecies);
    for (int i = 0; i < nspecies; ++i) {
        allcontainers[i] = std::make_unique<WarpXFluidContainer>(nlevs_max, i, species_names[i], geom);
    }
}

void
MultiFluidContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        const ParmParse pp_fluids("fluids");
        pp_fluids.queryarr("species_names", species_names);

        initialized = true;
    }
}

WarpXFluidContainer&
MultiFluidContainer::GetFluidContainerFromName (const std::string& name) const
{
    auto it = std::find(species_names.begin(), species_names.end(), name);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        it != species_names.end(),
        "unknown species name");
    const int i = std::distance(species_names.begin(), it);
    return *allcontainers[i];
}

void
MultiFluidContainer::AllocateLevelMFs (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    for (auto& fl : allcontainers) {
        fl->AllocateLevelMFs(lev, ba, dm);
    }
}

void
MultiFluidContainer::InitData (int lev, amrex::Box init_box, amrex::Real cur_time)
{
    for (auto& fl : allcontainers) {
        fl->InitData(lev, init_box, cur_time);
    }
}


void
MultiFluidContainer::DepositCharge (int lev, amrex::MultiFab &rho)
{
    for (auto& fl : allcontainers) {
        fl->DepositCharge(lev,rho);
    }
}

void
MultiFluidContainer::Evolve (int lev,
                            const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                            const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                            MultiFab* rho, MultiFab& jx, MultiFab& jy, MultiFab& jz,
                            amrex::Real cur_time, bool skip_deposition)
{
    for (auto& fl : allcontainers) {
        fl->Evolve(lev, Ex, Ey, Ez, Bx, By, Bz, rho, jx, jy, jz, cur_time, skip_deposition);
    }
}
