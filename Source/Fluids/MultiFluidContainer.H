/* Copyright 2023 Grant Johnson, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MultiFluidContainer_H_
#define WARPX_MultiFluidContainer_H_
#include "Evolve/WarpXDtType.H"

#include "WarpXFluidContainer_fwd.H"

#include<AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <string>

/**
 * The class MultiFluidContainer holds multiple instances of the
 * class WarpXFluidContainer, stored in its member variable "allcontainers".
 * The class WarpX typically has a single (pointer to an) instance of
 * MultiFluidContainer.
 *
 * MultiFluidContainer typically has two types of functions:
 * - Functions that loop over all instances of WarpXFluidContainer in
 *   allcontainers and calls the corresponding function (for instance,
 *   MultiFluidContainer::Evolve loops over all fluid containers and
 *   calls the corresponding WarpXFluidContainer::Evolve function).
 * - Functions that specifically handle multiple species (for instance
 *   ReadParameters).
 */
class MultiFluidContainer
{

public:

    MultiFluidContainer (int nlevs_max);

    ~MultiFluidContainer() = default;

    MultiFluidContainer (MultiFluidContainer const &)              = delete;
    MultiFluidContainer& operator= (MultiFluidContainer const & )  = delete;
    MultiFluidContainer(MultiFluidContainer&& )                    = default;
    MultiFluidContainer& operator=(MultiFluidContainer&& )         = default;

    [[nodiscard]] WarpXFluidContainer&
    GetFluidContainer (int ispecies) const {return *allcontainers[ispecies];}

#ifdef WARPX_USE_OPENPMD
    std::unique_ptr<WarpXFluidContainer>& GetUniqueContainer(int ispecies) {
      return  allcontainers[ispecies];
    }
#endif

    void AllocateLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void InitData (int lev, amrex::Box init_box, amrex::Real cur_time);

    ///
    /// This evolves all the fluids by one PIC time step, including current deposition, the
    /// field solve, and pushing the fluids, for all the species in the MultiFluidContainer.
    ///
    void Evolve (int lev,
                 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
                 amrex::MultiFab* rho, amrex::MultiFab& jx,  amrex::MultiFab& jy, amrex::MultiFab& jz,
                 amrex::Real cur_time, bool skip_deposition=false);

    [[nodiscard]] int nSpecies() const {return static_cast<int>(species_names.size());}

    void DepositCharge (int lev, amrex::MultiFab &rho);
    void DepositCurrent (int lev,
        amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz);

private:

    std::vector<std::string> species_names;

    // Vector of fluid species
    amrex::Vector<std::unique_ptr<WarpXFluidContainer>> allcontainers;

};
#endif /*WARPX_MultiFluidContainer_H_*/
