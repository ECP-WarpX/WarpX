/* Copyright 2023 Grant Johnson, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_WarpXFluidContainer_H_
#define WARPX_WarpXFluidContainer_H_

#include "Evolve/WarpXDtType.H"
#include "Initialization/PlasmaInjector.H"
#include "MultiFluidContainer.H"

#include<AMReX_MultiFab.H>
#include<AMReX_Vector.H>

#include <string>

using namespace amrex::literals;

/**
 * WarpXFluidContainer is the base polymorphic class from which all concrete
 * fluid container classes derive.
 *
 * WarpXFluidContainer contains the main functions for initialization,
 * interaction with the grid (field gather and current deposition), fluid
 * source and push.
 */
class WarpXFluidContainer
{
public:
    friend MultiFluidContainer;

    WarpXFluidContainer (int nlevs_max, int ispecies, const std::string& name);
    ~WarpXFluidContainer() {}

    void AllocateLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void InitData (int lev);

    void ReadParameters ();
    /**
     * Evolve updates a single timestep (dt) of the cold relativistic fluid eqautions
     */
    void Evolve (int lev,
        const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
        const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
        amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, bool skip_deposition=false);

    /**
     * AdvectivePush_Muscl takes a single timestep (dt) of the cold relativistic fluid equations
     * using a Muscl-Handcock scheme
     */
    /**
     * \brief Advective term, cold-rel. fluids
     *
     * \param[in] lev refinement level
     */
    void AdvectivePush_Muscl (int lev);

    /**
     * centrifugal_source adds contributions due to curvature acceleration for a
     * single timestep using an SSP-RK3 timestep
    */
    /**
     * \brief Centrifugal source term
     *
     * \param[in] lev refinement level
     */
    void centrifugal_source (int lev);

    /**
     * GatherAndPush introduces the lorentz term in the cold relativistic fluid
     * equations for a single timestep (dt) using Higuera and Cary Push
     */
    /**
     * \brief Lorentz Momentum Source
     *
     * \param[in] lev refinement level
     * \param[in] Ex Yee electric field (x)
     * \param[in] Ey Yee electric field (y)
     * \param[in] Ez Yee electric field (z)
     * \param[in] Bx Yee magnetic field (x)
     * \param[in] By Yee magnetic field (y)
     * \param[in] Bz Yee magnetic field (z)
     */
    void GatherAndPush (int lev,
        const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
        const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz);

    /**
     * DepositCurrent interpolates the fluid current density comps. onto the Yee grid and
     * sums the contributions to the particle current density
     */
    /**
     * \brief Deposit fluid current density.
     *
     * \param[in] lev refinement level
     * \param[in,out] jx current density MultiFab x comp.
     * \param[in,out] jy current density MultiFab y comp.
     * \param[in,out] jz current density MultiFab z comp.
     */
    void DepositCurrent (int lev,
        amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz);

    /**
     * DepositCharge deposits density onto rho
     */
    /**
     * \brief Deposit fluid charge density.
     *
     * \param[in] lev refinement level
     * \param[in,out] rho charge density MultiFab.
     */
    void DepositCharge (int lev, amrex::MultiFab &rho);

    //amrex::Real getCharge () {return charge;}
    amrex::Real getCharge () const {return charge;}
    //amrex::Real getMass () {return mass;}
    amrex::Real getMass () const {return mass;}

protected:
    int species_id;
    std::string species_name;
    amrex::Real charge;
    amrex::Real mass;

    int do_not_push = 0;
    int do_not_gather = 0;
    int do_not_deposit = 0;
    std::string physical_element;
    PhysicalSpecies physical_species;

    std::unique_ptr<PlasmaInjector> plasma_injector;

private:
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > N;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > NU;
};

#endif
