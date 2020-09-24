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
#include "BoundaryConditions/PML.H"
#include "BoundaryConditions/PML_current.H"
#include "BoundaryConditions/PMLComponent.H"
#include "Utils/WarpXConst.H"
#include <AMReX_Gpu.H>
#include <AMReX.H>

using namespace amrex;

/**
 * \brief Update the E field, over one timestep
 */
void FiniteDifferenceSolver::EvolveEPML (
    std::array< amrex::MultiFab*, 3 > Efield,
    std::array< amrex::MultiFab*, 3 > const Bfield,
    std::array< amrex::MultiFab*, 3 > const Jfield,
    amrex::MultiFab* const Ffield,
    MultiSigmaBox const& sigba,
    amrex::Real const dt, bool pml_has_particles ) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    amrex::ignore_unused(Efield, Bfield, Jfield, Ffield, sigba, dt, pml_has_particles);
    amrex::Abort("PML are not implemented in cylindrical geometry.");
#else
    if (m_do_nodal) {

        EvolveEPMLCartesian <CartesianNodalAlgorithm> (
            Efield, Bfield, Jfield, Ffield, sigba, dt, pml_has_particles );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {

        EvolveEPMLCartesian <CartesianYeeAlgorithm> (
            Efield, Bfield, Jfield, Ffield, sigba, dt, pml_has_particles );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        EvolveEPMLCartesian <CartesianCKCAlgorithm> (
            Efield, Bfield, Jfield, Ffield, sigba, dt, pml_has_particles );

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
    MultiSigmaBox const& sigba,
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
        Box const& tex  = mfi.tilebox(Efield[0]->ixType().ixType());
        Box const& tey  = mfi.tilebox(Efield[1]->ixType().ixType());
        Box const& tez  = mfi.tilebox(Efield[2]->ixType().ixType());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tex, tey, tez,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ex(i, j, k, PMLComp::xz) -= c2 * dt * (
                    T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yx)
                  + T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yz) );
                Ex(i, j, k, PMLComp::xy) += c2 * dt * (
                    T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zx)
                  + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zy) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ey(i, j, k, PMLComp::yx) -= c2 * dt * (
                    T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k, PMLComp::zx)
                  + T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k, PMLComp::zy) );
                Ey(i, j, k, PMLComp::yz) += c2 * dt * (
                    T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k, PMLComp::xy)
                  + T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k, PMLComp::xz) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ez(i, j, k, PMLComp::zy) -= c2 * dt * (
                    T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k, PMLComp::xy)
                  + T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k, PMLComp::xz) );
                Ez(i, j, k, PMLComp::zx) += c2 * dt * (
                    T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k, PMLComp::yx)
                  + T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k, PMLComp::yz) );
            }

        );

        // If F is not a null pointer, further update E using the grad(F) term
        // (hyperbolic correction for errors in charge conservation)
        if (Ffield) {
            // Extract field data for this grid/tile
            Array4<Real> const& F = Ffield->array(mfi);

            // Loop over the cells and update the fields
            amrex::ParallelFor(tex, tey, tez,

                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Ex(i, j, k, PMLComp::xx) += c2 * dt * (
                        T_Algo::UpwardDx(F, coefs_x, n_coefs_x, i, j, k, PMLComp::x)
                      + T_Algo::UpwardDx(F, coefs_x, n_coefs_x, i, j, k, PMLComp::y)
                      + T_Algo::UpwardDx(F, coefs_x, n_coefs_x, i, j, k, PMLComp::z) );
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Ey(i, j, k, PMLComp::yy) += c2 * dt * (
                        T_Algo::UpwardDy(F, coefs_y, n_coefs_y, i, j, k, PMLComp::x)
                      + T_Algo::UpwardDy(F, coefs_y, n_coefs_y, i, j, k, PMLComp::y)
                      + T_Algo::UpwardDy(F, coefs_y, n_coefs_y, i, j, k, PMLComp::z) );
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k){
                    Ez(i, j, k, PMLComp::zz) += c2 * dt * (
                        T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, k, PMLComp::x)
                      + T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, k, PMLComp::y)
                      + T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, k, PMLComp::z) );
                }
            );
        }

        // Update the E field in the PML, using the current
        // deposited by the particles in the PML
        if (pml_has_particles) {

            // Extract field data for this grid/tile
            Array4<Real> const& Jx = Jfield[0]->array(mfi);
            Array4<Real> const& Jy = Jfield[1]->array(mfi);
            Array4<Real> const& Jz = Jfield[2]->array(mfi);
            const Real* sigmaj_x = sigba[mfi].sigma[0].data();
            const Real* sigmaj_y = sigba[mfi].sigma[1].data();
            const Real* sigmaj_z = sigba[mfi].sigma[2].data();
            int const x_lo = sigba[mfi].sigma[0].lo();
#if (AMREX_SPACEDIM == 3)
            int const y_lo = sigba[mfi].sigma[1].lo();
            int const z_lo = sigba[mfi].sigma[2].lo();
#else
            int const y_lo = 0;
            int const z_lo = sigba[mfi].sigma[1].lo();
#endif
            const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;

            amrex::ParallelFor( tex, tey, tez,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    push_ex_pml_current(i, j, k, Ex, Jx,
                        sigmaj_y, sigmaj_z, y_lo, z_lo, mu_c2_dt);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    push_ey_pml_current(i, j, k, Ey, Jy,
                        sigmaj_x, sigmaj_z, x_lo, z_lo, mu_c2_dt);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    push_ez_pml_current(i, j, k, Ez, Jz,
                        sigmaj_x, sigmaj_y, x_lo, y_lo, mu_c2_dt);
                }
            );
        }

    }

}

#endif // corresponds to ifndef WARPX_DIM_RZ
