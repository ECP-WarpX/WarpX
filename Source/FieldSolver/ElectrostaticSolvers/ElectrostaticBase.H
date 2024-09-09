/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenwald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */

#include "PoissonBoundaryHandler.H"
#include "WarpX.H"

#include <AMReX_Array.H>

using namespace amrex;

class ElectrostaticBase
{
public:
    void ElectrostaticBase( int nlevs_max );

    virtual void AllocateMFs ( int nlevs_max ) {};

    void setPhiBC (amrex::Vector<std::unique_ptr<amrex::MultiFab>>& phi ) const;

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

    PoissonBoundaryHandler m_poisson_boundary_handler;
}
