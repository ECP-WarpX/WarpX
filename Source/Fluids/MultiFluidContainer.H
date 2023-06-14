/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Jean-Luc Vay, Junmin Gu
 * Luca Fedeli, Mathieu Lobet, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MultiFluidContainer_H_
#define WARPX_MultiFluidContainer_H_

#include "WarpXFluidContainer_fwd.H"

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_GpuControl.H>
#include <AMReX_INT.H>
#include <AMReX_MFIter.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>
#include <AMReX_AmrCoreFwd.H>

#include <algorithm>
#include <array>
#include <iosfwd>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <vector>

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

    MultiFluidContainer (amrex::AmrCore* amr_core);

    ~MultiFluidContainer() {}

    WarpXFluidContainer&
    GetFluidContainer (int ispecies) const {return *allcontainers[ispecies];}

    WarpXFluidContainer*
    GetFluidContainerPtr (int ispecies) const {return allcontainers[ispecies].get();}

    WarpXFluidContainer&
    GetFluidContainerFromName (const std::string& name) const;

#ifdef WARPX_USE_OPENPMD
    std::unique_ptr<WarpXFluidContainer>& GetUniqueContainer(int ispecies) {
      return  allcontainers[ispecies];
    }
#endif

    void AllocData ();

    void InitData ();

    ///
    /// This evolves all the fluids by one PIC time step, including current deposition, the
    /// field solve, and pushing the fluids, for all the species in the MultiFluidContainer.
    ///
    /*  void Evolve (int lev,
                 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
                 amrex::MultiFab& jx,  amrex::MultiFab& jy, amrex::MultiFab& jz,
                 amrex::MultiFab* cjx,  amrex::MultiFab* cjy, amrex::MultiFab* cjz,
                 amrex::MultiFab* rho, amrex::MultiFab* crho,
                 const amrex::MultiFab* cEx, const amrex::MultiFab* cEy, const amrex::MultiFab* cEz,
                 const amrex::MultiFab* cBx, const amrex::MultiFab* cBy, const amrex::MultiFab* cBz,
                 amrex::Real t, amrex::Real dt, DtType a_dt_type=DtType::Full, bool skip_deposition=false); */

    int nSpecies() const {return species_names.size();}

    std::vector<std::string> GetSpeciesNames() const { return species_names; }

protected:

    std::vector<std::string> species_names;

private:

    // Vector of fluid species
    amrex::Vector<std::unique_ptr<WarpXFluidContainer>> allcontainers;

    void ReadParameters ();

};
#endif /*WARPX_MultiFluidContainer_H_*/
