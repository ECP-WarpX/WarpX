/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Utils/WarpXAlgorithmSelection.H"
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#ifdef WARPX_DIM_RZ
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#endif
#include "BoundaryConditions/PMLComponent.H"
#include "Utils/WarpXConst.H"
#include <AMReX_Gpu.H>

using namespace amrex;

/**
 * \brief Update the E field, over one timestep
 */
void FiniteDifferenceSolver::EvolveEPML (
    std::array< amrex::MultiFab*, 3 > Efield,
    std::array< amrex::MultiFab*, 3 > const Bfield,
    std::array< amrex::MultiFab*, 3 > const Jfield,
    amrex::MultiFab* const Ffield,
    amrex::Real const dt, bool pml_has_particles ) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    amrex::Abort("PML are not implemented in cylindrical geometry.");
#else
    if (m_do_nodal) {

        EvolveEPMLCartesian <CartesianNodalAlgorithm> (
            Efield, Bfield, Jfield, Ffield, dt, pml_has_particles );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {

        EvolveEPMLCartesian <CartesianYeeAlgorithm> (
            Efield, Bfield, Jfield, Ffield, dt, pml_has_particles );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        EvolveEPMLCartesian <CartesianCKCAlgorithm> (
            Efield, Bfield, Jfield, Ffield, dt, pml_has_particles );

    } else {
        amrex::Abort("Unknown algorithm");
    }
#endif
}


#ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveEPMLCartesian (
    std::array< amrex::MultiFab*, 3 > Efield,
    std::array< amrex::MultiFab*, 3 > const Bfield,
    std::array< amrex::MultiFab*, 3 > const Jfield,
    amrex::MultiFab* const Ffield,
    amrex::Real const dt, bool pml_has_particles ) {

    Real constexpr c2 = PhysConst::c * PhysConst::c;

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        Array4<Real> const& Ex = Efield[0]->array(mfi);
        Array4<Real> const& Ey = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);
        Array4<Real> const& Bx = Bfield[0]->array(mfi);
        Array4<Real> const& By = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        int const n_coefs_x = m_stencil_coefs_x.size();
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        int const n_coefs_y = m_stencil_coefs_y.size();
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        int const n_coefs_z = m_stencil_coefs_z.size();

        // Extract tileboxes for which to loop
        Box const& tbx  = mfi.tilebox(Efield[0]->ixType().ixType());
        Box const& tby  = mfi.tilebox(Efield[1]->ixType().ixType());
        Box const& tbz  = mfi.tilebox(Efield[2]->ixType().ixType());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tbx, tby, tbz,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ex(i, j, k, PMLComp::xz) -= c2 * dt * (
                    T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yx)
                  + T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yy)
                  + T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yz) );
                Ex(i, j, k, PMLComp::xy) += c2 * dt * (
                    T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zx)
                  + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zy)
                  + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zz) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ey(i, j, k, PMLComp::yx) -= c2 * dt * (
                    T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k, PMLComp::zx)
                  + T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k, PMLComp::zy)
                  + T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k, PMLComp::zz) );
                Ey(i, j, k, PMLComp::yz) += c2 * dt * (
                    T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k, PMLComp::xx)
                  + T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k, PMLComp::xy)
                  + T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k, PMLComp::xz) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ez(i, j, k, PMLComp::zy) -= c2 * dt * (
                    T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k, PMLComp::xx)
                  + T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k, PMLComp::xy)
                  + T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k, PMLComp::xz) );
                Ez(i, j, k, PMLComp::zx) += c2 * dt * (
                    T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k, PMLComp::yx)
                  + T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k, PMLComp::yy)
                  + T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k, PMLComp::yz) );
            }

        );

    }

}

#endif // corresponds to ifndef WARPX_DIM_RZ
