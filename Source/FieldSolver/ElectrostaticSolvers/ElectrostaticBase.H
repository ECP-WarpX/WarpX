/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ELECTROSTATIC_BASE_H_
#define WARPX_ELECTROSTATIC_BASE_H_

#include "PoissonBoundaryHandler.H"
#include "Fluids/MultiFluidContainer.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX_Array.H>

using namespace amrex;

class ElectrostaticBase
{
public:
    ElectrostaticBase() = default;
    ElectrostaticBase( int nlevs_max );

    virtual ~ElectrostaticBase();

    // Prohibit Move and Copy operations
    ElectrostaticBase(const ElectrostaticBase&) = delete;
    ElectrostaticBase& operator=(const ElectrostaticBase&) = delete;
    ElectrostaticBase(ElectrostaticBase&&) = delete;
    ElectrostaticBase& operator=(ElectrostaticBase&&) = delete;


    virtual void ReadParameters () {};
    virtual void AllocateMFs ( int nlevs_max ) {};

    void setPhiBC (
        amrex::Vector<std::unique_ptr<amrex::MultiFab>>& phi,
        amrex::Real t
    ) const;

    void computePhi (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
        std::array<Real, 3> const beta,
        Real const required_precision,
        Real absolute_tolerance,
        int const max_iters,
        int const verbosity
    ) const;

    void computeE (
        amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
        std::array<amrex::Real, 3> const beta
    ) const;
    void computeB (
        amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& B,
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
        std::array<amrex::Real, 3> const beta
    ) const;

    int max_level;
    std::unique_ptr<PoissonBoundaryHandler> m_poisson_boundary_handler;

    // Parameters for lab frame electrostatic
    amrex::Real self_fields_required_precision = 1e-11_rt;
    amrex::Real self_fields_absolute_tolerance = 0.0_rt;
    int self_fields_max_iters = 200;
    int self_fields_verbosity = 2;
};

#endif
