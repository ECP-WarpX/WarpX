/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote
 * Jean-Luc Vay, Junmin Gu, Luca Fedeli
 * Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 * Weiqun Zhang, Yinjian Zhao
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
#include <AMReX_Vector.H>

#include <string>

using namespace amrex::literals;

/**
 * TODO
 */
class WarpXFluidContainer
{
public:
    friend MultiFluidContainer;

    WarpXFluidContainer (int nlevs_max, int ispecies, const std::string& name);
    ~WarpXFluidContainer() {}

    void AllocateLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    void InitData () {}

    void ReadParameters ();
    /**
     * Evolve is the central WarpXFluidContainer function that advances
     * fluids for a time dt (typically one timestep).
     */
    void Evolve (int lev,
        const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
        const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
        amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz,
        amrex::MultiFab* cjx, amrex::MultiFab* cjy, amrex::MultiFab* cjz,
        amrex::MultiFab* rho, amrex::MultiFab* crho,
        const amrex::MultiFab* cEx, const amrex::MultiFab* cEy, const amrex::MultiFab* cEz,
        const amrex::MultiFab* cBx, const amrex::MultiFab* cBy, const amrex::MultiFab* cBz,
        amrex::Real t, amrex::Real dt, DtType a_dt_type=DtType::Full, bool skip_deposition=false) {}

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
