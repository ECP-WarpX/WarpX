/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"

#include "BoundaryConditions/PML.H"
#include "BoundaryConditions/PMLComponent.H"
#include "BoundaryConditions/PML_current.H"

#ifndef WARPX_DIM_RZ
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#else
#   include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#endif
#include "Utils/TextMsg.H"
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

using namespace amrex;

/**
 * \brief Update the B field, over one timestep
 */

// Taking the simplest case (v = c) and the angle here and in the input file for first test
Real constexpr c2 = PhysConst::c * PhysConst::c;
Real constexpr c2v = PhysConst::c;
Real thet_pml = MathConst::pi/6;
Real tn = std::tan(thet_pml);
Real sn = std::sin(thet_pml);
Real cs = std::cos(thet_pml);

void FiniteDifferenceSolver::EvolveBPML (
    std::array< amrex::MultiFab*, 3 > Bfield,
    std::array< amrex::MultiFab*, 3 > const Efield,
    amrex::Real const dt,
    const bool dive_cleaning,
    std::array< amrex::MultiFab*, 3 > const Jfield,
    MultiSigmaBox const& sigba, bool pml_has_particles) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    amrex::ignore_unused(Bfield, Efield, dt, dive_cleaning);
    amrex::Abort(Utils::TextMsg::Err(
        "PML are not implemented in cylindrical geometry."));
#else
    if (m_do_nodal) {

        EvolveBPMLCartesian <CartesianNodalAlgorithm> (Bfield, Efield, dt, dive_cleaning, Jfield, sigba, pml_has_particles);

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee || m_fdtd_algo == MaxwellSolverAlgo::ECT) {

        EvolveBPMLCartesian <CartesianYeeAlgorithm> (Bfield, Efield, dt, dive_cleaning, Jfield, sigba, pml_has_particles);

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        EvolveBPMLCartesian <CartesianCKCAlgorithm> (Bfield, Efield, dt, dive_cleaning, Jfield, sigba, pml_has_particles);

    } else {
        amrex::Abort("EvolveBPML: Unknown algorithm");
    }
#endif
}


#ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveBPMLCartesian (
    std::array< amrex::MultiFab*, 3 > Bfield,
    std::array< amrex::MultiFab*, 3 > const Efield,
    amrex::Real const dt,
    const bool dive_cleaning,
    std::array< amrex::MultiFab*, 3 > const Jfield,
    MultiSigmaBox const& sigba, bool pml_has_particles) {

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bfield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        Array4<Real> const& Bx = Bfield[0]->array(mfi);
        Array4<Real> const& By = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);
        Array4<Real> const& Ex = Efield[0]->array(mfi);
        Array4<Real> const& Ey = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        int const n_coefs_x = m_stencil_coefs_x.size();
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        int const n_coefs_y = m_stencil_coefs_y.size();
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        int const n_coefs_z = m_stencil_coefs_z.size();

        // My addition
        const Real* sigma_x = sigba[mfi].sigma[0].data();
        const Real* sigma_y = sigba[mfi].sigma[1].data();
        const Real* sigma_z = sigba[mfi].sigma[2].data();

        int const x_lo = sigba[mfi].sigma[0].lo();
#if defined(WARPX_DIM_3D)
        int const y_lo = sigba[mfi].sigma[1].lo();
        int const z_lo = sigba[mfi].sigma[2].lo();
#else
        int const y_lo = 0;
        int const z_lo = sigba[mfi].sigma[1].lo();
#endif
        const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;

        // Extract tileboxes for which to loop
        Box const& tbx  = mfi.tilebox(Bfield[0]->ixType().ixType());
        Box const& tby  = mfi.tilebox(Bfield[1]->ixType().ixType());
        Box const& tbz  = mfi.tilebox(Bfield[2]->ixType().ixType());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tbx, tby, tbz,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                amrex::Real UpwardDz_Ey_yy = 0._rt;
                amrex::Real UpwardDy_Ez_zz = 0._rt;
                amrex::Real UpwardDy_By_yy = 0._rt;
                if (dive_cleaning)
                {
                    UpwardDz_Ey_yy = T_Algo::UpwardDz(Ey, coefs_z, n_coefs_z, i, j, k, PMLComp::yy);
                    UpwardDy_Ez_zz = T_Algo::UpwardDy(Ez, coefs_y, n_coefs_y, i, j, k, PMLComp::zz);
                    UpwardDy_By_yy = T_Algo::UpwardDy(By, coefs_y, n_coefs_y, i, j, k, PMLComp::yy);
                }

                Bx(i, j, k, PMLComp::xz) += dt * (
                    T_Algo::UpwardDz(Ey, coefs_z, n_coefs_z, i, j, k, PMLComp::yx)
                  + T_Algo::UpwardDz(Ey, coefs_z, n_coefs_z, i, j, k, PMLComp::yz)
                  + UpwardDz_Ey_yy
                  + c2v*sn*(T_Algo::UpwardDy(By, coefs_y, n_coefs_y, i, j, k, PMLComp::yz)
                  + T_Algo::UpwardDy(By, coefs_y, n_coefs_y, i, j, k, PMLComp::yx)
                  + UpwardDy_By_yy));

                Bx(i, j, k, PMLComp::xy) -= dt * (
                    T_Algo::UpwardDy(Ez, coefs_y, n_coefs_y, i, j, k, PMLComp::zx)
                  + T_Algo::UpwardDy(Ez, coefs_y, n_coefs_y, i, j, k, PMLComp::zy)
                  + UpwardDy_Ez_zz
                  + c2v*sn*(T_Algo::UpwardDy(By, coefs_y, n_coefs_y, i, j, k, PMLComp::yz)
                  + T_Algo::UpwardDy(By, coefs_y, n_coefs_y, i, j, k, PMLComp::yx)
                  + UpwardDy_By_yy));
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                amrex::Real UpwardDx_Ez_zz = 0._rt;
                amrex::Real UpwardDz_Ex_xx = 0._rt;
                amrex::Real UpwardDy_Bz_zz = 0._rt;
                amrex::Real UpwardDx_Ex_xx = 0._rt;
                amrex::Real UpwardDz_By_yy = 0._rt;

                const Real sigma_xy = sigma_y[j-y_lo];
                const Real sigma_xz = sigma_z[k-z_lo];

                if (dive_cleaning)
                {
                    UpwardDx_Ez_zz = T_Algo::UpwardDx(Ez, coefs_x, n_coefs_x, i, j, k, PMLComp::zz);
                    UpwardDz_Ex_xx = T_Algo::UpwardDz(Ex, coefs_z, n_coefs_z, i, j, k, PMLComp::xx);
                    UpwardDy_Bz_zz = T_Algo::UpwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zz);
                    UpwardDx_Ex_xx = T_Algo::UpwardDz(Ex, coefs_x, n_coefs_x, i, j, k, PMLComp::xx);
                    UpwardDz_By_yy = T_Algo::UpwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yy);
                }

                By(i, j, k, PMLComp::yx) += dt * (
                    T_Algo::UpwardDx(Ez, coefs_x, n_coefs_x, i, j, k, PMLComp::zx)
                  + T_Algo::UpwardDx(Ez, coefs_x, n_coefs_x, i, j, k, PMLComp::zy)
                  + UpwardDx_Ez_zz
                  + tn*(T_Algo::UpwardDx(Ex, coefs_x, n_coefs_x, i, j, k, PMLComp::xz)
                  + T_Algo::UpwardDx(Ex, coefs_x, n_coefs_x, i, j, k, PMLComp::xy)
                  + UpwardDx_Ex_xx)
                  + c2v*tn*sn*(T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zx)
                  + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zy)
                  + UpwardDy_Bz_zz
                  - T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yz)
                  - T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yx)
                  - UpwardDz_By_yy)
                  - (tn*sn/c2v)*sigma_xz*mu_c2_dt*Ex(i, j, k, PMLComp::xz));

                By(i, j, k, PMLComp::yz) -= dt * (
                    UpwardDz_Ex_xx
                  + T_Algo::UpwardDz(Ex, coefs_z, n_coefs_z, i, j, k, PMLComp::xy)
                  + T_Algo::UpwardDz(Ex, coefs_z, n_coefs_z, i, j, k, PMLComp::xz)
                  - c2v*tn*sn*(UpwardDz_By_yy
                  + T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yz)
                  + T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k, PMLComp::yx))
                  + tn*(T_Algo::UpwardDx(Ex, coefs_x, n_coefs_x, i, j, k, PMLComp::xz)
                  + T_Algo::UpwardDx(Ex, coefs_x, n_coefs_x, i, j, k, PMLComp::xy)
                  + UpwardDx_Ex_xx)
                  + c2v*tn*sn*(T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zx)
                  + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k, PMLComp::zy)
                  + UpwardDy_Bz_zz)
                  - (tn*sn/c2v)*sigma_xz*mu_c2_dt*Ex(i, j, k, PMLComp::xz));
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                amrex::Real UpwardDy_Ex_xx = 0._rt;
                amrex::Real UpwardDx_Ey_yy = 0._rt;
                if (dive_cleaning)
                {
                    UpwardDy_Ex_xx = T_Algo::UpwardDy(Ex, coefs_y, n_coefs_y, i, j, k, PMLComp::xx);
                    UpwardDx_Ey_yy = T_Algo::UpwardDx(Ey, coefs_x, n_coefs_x, i, j, k, PMLComp::yy);
                }

                Bz(i, j, k, PMLComp::zy) += dt * (
                    UpwardDy_Ex_xx
                  + T_Algo::UpwardDy(Ex, coefs_y, n_coefs_y, i, j, k, PMLComp::xy)
                  + T_Algo::UpwardDy(Ex, coefs_y, n_coefs_y, i, j, k, PMLComp::xz) );

                Bz(i, j, k, PMLComp::zx) -= dt * (
                    T_Algo::UpwardDx(Ey, coefs_x, n_coefs_x, i, j, k, PMLComp::yx)
                  + T_Algo::UpwardDx(Ey, coefs_x, n_coefs_x, i, j, k, PMLComp::yz)
                  + UpwardDx_Ey_yy);
            }

        );

        // Update the B field in the PML, using the current
        // deposited by the particles in the PML

        if (pml_has_particles) {

            // Extract field data for this grid/tile
            Array4<Real> const& Jx = Jfield[0]->array(mfi);
            const Real* sigmaj_x = sigba[mfi].sigma[0].data();
            const Real* sigmaj_y = sigba[mfi].sigma[1].data();
            const Real* sigmaj_z = sigba[mfi].sigma[2].data();
            int const x_lo = sigba[mfi].sigma[0].lo();
#if defined(WARPX_DIM_3D)
            int const y_lo = sigba[mfi].sigma[1].lo();
            int const z_lo = sigba[mfi].sigma[2].lo();
#else
            int const y_lo = 0;
            int const z_lo = sigba[mfi].sigma[1].lo();
#endif
            const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;
            const Real C_th = tn*sn/c2v;

        amrex::ParallelFor(tby,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    push_by_pml_current(i, j, k, By, Jx,
                                sigmaj_y, sigmaj_z, y_lo, z_lo, mu_c2_dt, C_th);
                }
            );

        }
    }
}

#endif // corresponds to ifndef WARPX_DIM_RZ
