/* Copyright 2020
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */


#include "FiniteDifferenceSolver.H"
#ifndef WARPX_DIM_RZ
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#else
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#endif
#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

using namespace amrex;

void FiniteDifferenceSolver::EvolveG (
    std::unique_ptr<amrex::MultiFab>& Gfield,
    std::array<std::unique_ptr<amrex::MultiFab>,3> const& Bfield,
    amrex::Real const dt)
{
#ifdef WARPX_DIM_RZ
    // TODO Implement G update equation in RZ geometry
    amrex::ignore_unused(Gfield, Bfield, dt);
#else
    // Select algorithm
    if (m_do_nodal)
    {
        EvolveGCartesian<CartesianNodalAlgorithm>(Gfield, Bfield, dt);
    }
    else if (m_fdtd_algo == MaxwellSolverAlgo::Yee)
    {
        EvolveGCartesian<CartesianYeeAlgorithm>(Gfield, Bfield, dt);
    }
    else if (m_fdtd_algo == MaxwellSolverAlgo::CKC)
    {
        EvolveGCartesian<CartesianCKCAlgorithm>(Gfield, Bfield, dt);
    }
    else
    {
        amrex::Abort("EvolveG: unknown FDTD algorithm");
    }
#endif
}

#ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveGCartesian (
    std::unique_ptr<amrex::MultiFab>& Gfield,
    std::array<std::unique_ptr<amrex::MultiFab>,3> const& Bfield,
    amrex::Real const dt)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

    // Loop over grids and over tiles within each grid
    for (amrex::MFIter mfi(*Gfield, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Extract field data for this grid/tile
        amrex::Array4<amrex::Real> const& G = Gfield->array(mfi);
        amrex::Array4<amrex::Real> const& Bx = Bfield[0]->array(mfi);
        amrex::Array4<amrex::Real> const& By = Bfield[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Bz = Bfield[2]->array(mfi);

        // Extract stencil coefficients
        amrex::Real const* const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        amrex::Real const* const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        amrex::Real const* const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();

        const int n_coefs_x = m_stencil_coefs_x.size();
        const int n_coefs_y = m_stencil_coefs_y.size();
        const int n_coefs_z = m_stencil_coefs_z.size();

        // Extract tilebox to loop over
        amrex::Box const& tf = mfi.tilebox(Gfield->ixType().toIntVect());

        // Loop over cells and update G
        amrex::ParallelFor(tf, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            G(i,j,k) += dt * (T_Algo::UpwardDx(Bx, coefs_x, n_coefs_x, i, j, k)
                            + T_Algo::UpwardDy(By, coefs_y, n_coefs_y, i, j, k)
                            + T_Algo::UpwardDz(Bz, coefs_z, n_coefs_z, i, j, k));
        });
    }
}

#endif
