/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Jean-Luc Vay, Luca Fedeli
 * Mathieu Lobet, Maxence Thevenet, Neil Zaim
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include "WarpX.H"

#include <ablastr/utils/Communication.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace amrex;

MultiFluidContainer::MultiFluidContainer (int nlevs_max)
{
    ReadParameters();

    auto const nspecies = static_cast<int>(species_names.size());

    allcontainers.resize(nspecies);
    for (int i = 0; i < nspecies; ++i) {
        allcontainers[i] = std::make_unique<WarpXFluidContainer>(nlevs_max, i, species_names[i]);
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
    for (auto& pc : allcontainers) {
        pc->AllocateLevelMFs(lev, ba, dm);
    }
}

void
MultiFluidContainer::InitData ()
{
    for (auto& pc : allcontainers) {
        pc->InitData();
    }
}


void
MultiFluidContainer::Evolve (int lev,
                            const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                            const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                            MultiFab& jx, MultiFab& jy, MultiFab& jz,
                            MultiFab* cjx,  MultiFab* cjy, MultiFab* cjz,
                            MultiFab* rho, MultiFab* crho,
                            const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                            const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                            Real t, Real dt, DtType a_dt_type, bool skip_deposition)
{
    for (auto& pc : allcontainers) {
        pc->Evolve(lev, Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, cjx, cjy, cjz,
                   rho, crho, cEx, cEy, cEz, cBx, cBy, cBz, t, dt, a_dt_type, skip_deposition);
    }
}
