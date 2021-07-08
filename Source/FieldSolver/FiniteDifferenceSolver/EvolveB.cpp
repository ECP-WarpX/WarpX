/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FiniteDifferenceSolver.H"

#ifndef WARPX_DIM_RZ
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#else
#   include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#endif
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

using namespace amrex;

/**
 * \brief Update the B field, over one timestep
 */
void FiniteDifferenceSolver::EvolveB (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
    std::unique_ptr<amrex::MultiFab> const& Gfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& face_areas,
    int lev, amrex::Real const dt ) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    if (m_fdtd_algo == MaxwellSolverAlgo::Yee){
        ignore_unused(Gfield, face_areas);
        EvolveBCylindrical <CylindricalYeeAlgorithm> ( Bfield, Efield, lev, dt );
#else
    if (m_do_nodal) {

        EvolveBCartesian <CartesianNodalAlgorithm> ( Bfield, Efield, Gfield, face_areas, lev, dt );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {

        EvolveBCartesian <CartesianYeeAlgorithm> ( Bfield, Efield, Gfield, face_areas, lev, dt );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        EvolveBCartesian <CartesianCKCAlgorithm> ( Bfield, Efield, Gfield, face_areas, lev, dt );

#endif
    } else {
        amrex::Abort("EvolveB: Unknown algorithm");
    }
}


#ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveBCartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
    std::unique_ptr<amrex::MultiFab> const& Gfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& face_areas,
    int lev, amrex::Real const dt ) {

#ifndef AMREX_USE_EB
    amrex::ignore_unused(face_areas);
#endif

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    Real constexpr c2 = PhysConst::c * PhysConst::c;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bfield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        // Extract field data for this grid/tile
        Array4<Real> const& Bx = Bfield[0]->array(mfi);
        Array4<Real> const& By = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);
        Array4<Real> const& Ex = Efield[0]->array(mfi);
        Array4<Real> const& Ey = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);

#ifdef AMREX_USE_EB
        amrex::Array4<amrex::Real> const& Sx = face_areas[0]->array(mfi);
        amrex::Array4<amrex::Real> const& Sy = face_areas[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Sz = face_areas[2]->array(mfi);
#endif

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        int const n_coefs_x = m_stencil_coefs_x.size();
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        int const n_coefs_y = m_stencil_coefs_y.size();
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        int const n_coefs_z = m_stencil_coefs_z.size();

        // Extract tileboxes for which to loop
        Box const& tbx  = mfi.tilebox(Bfield[0]->ixType().toIntVect());
        Box const& tby  = mfi.tilebox(Bfield[1]->ixType().toIntVect());
        Box const& tbz  = mfi.tilebox(Bfield[2]->ixType().toIntVect());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tbx, tby, tbz,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip field push if this cell is fully covered by embedded boundaries
                if (Sx(i, j, k) <= 0) return;
#endif
                Bx(i, j, k) += dt * T_Algo::UpwardDz(Ey, coefs_z, n_coefs_z, i, j, k)
                             - dt * T_Algo::UpwardDy(Ez, coefs_y, n_coefs_y, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip field push if this cell is fully covered by embedded boundaries
                if (Sy(i, j, k) <= 0) return;
#endif
                By(i, j, k) += dt * T_Algo::UpwardDx(Ez, coefs_x, n_coefs_x, i, j, k)
                             - dt * T_Algo::UpwardDz(Ex, coefs_z, n_coefs_z, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip field push if this cell is fully covered by embedded boundaries
                if (Sz(i, j, k) <= 0) return;
#endif
                Bz(i, j, k) += dt * T_Algo::UpwardDy(Ex, coefs_y, n_coefs_y, i, j, k)
                             - dt * T_Algo::UpwardDx(Ey, coefs_x, n_coefs_x, i, j, k);
            }
        );

        // div(B) cleaning correction for errors in magnetic Gauss law (div(B) = 0)
        if (Gfield)
        {
            // Extract field data for this grid/tile
            Array4<Real> G = Gfield->array(mfi);

            // Loop over cells and update G
            amrex::ParallelFor(tbx, tby, tbz,

                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Bx(i,j,k) += c2 * dt * T_Algo::DownwardDx(G, coefs_x, n_coefs_x, i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    By(i,j,k) += c2 * dt * T_Algo::DownwardDy(G, coefs_y, n_coefs_y, i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Bz(i,j,k) += c2 * dt * T_Algo::DownwardDz(G, coefs_z, n_coefs_z, i, j, k);
                }
            );
        }

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }
}

#else // corresponds to ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveBCylindrical (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
    int lev, amrex::Real const dt ) {

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bfield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        // Extract field data for this grid/tile
        Array4<Real> const& Br = Bfield[0]->array(mfi);
        Array4<Real> const& Bt = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);
        Array4<Real> const& Er = Efield[0]->array(mfi);
        Array4<Real> const& Et = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_r = m_stencil_coefs_r.dataPtr();
        int const n_coefs_r = m_stencil_coefs_r.size();
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        int const n_coefs_z = m_stencil_coefs_z.size();

        // Extract cylindrical specific parameters
        Real const dr = m_dr;
        int const nmodes = m_nmodes;
        Real const rmin = m_rmin;

        // Extract tileboxes for which to loop
        Box const& tbr  = mfi.tilebox(Bfield[0]->ixType().toIntVect());
        Box const& tbt  = mfi.tilebox(Bfield[1]->ixType().toIntVect());
        Box const& tbz  = mfi.tilebox(Bfield[2]->ixType().toIntVect());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tbr, tbt, tbz,

            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){
                Real const r = rmin + i*dr; // r on nodal point (Br is nodal in r)
                if (r != 0) { // Off-axis, regular Maxwell equations
                    Br(i, j, 0, 0) += dt * T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 0); // Mode m=0
                    for (int m=1; m<nmodes; m++) { // Higher-order modes
                        Br(i, j, 0, 2*m-1) += dt*(
                            T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m-1)
                            - m * Ez(i, j, 0, 2*m  )/r );  // Real part
                        Br(i, j, 0, 2*m  ) += dt*(
                            T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m  )
                            + m * Ez(i, j, 0, 2*m-1)/r ); // Imaginary part
                    }
                } else { // r==0: On-axis corrections
                    // Ensure that Br remains 0 on axis (except for m=1)
                    Br(i, j, 0, 0) = 0.; // Mode m=0
                    for (int m=1; m<nmodes; m++) { // Higher-order modes
                        if (m == 1){
                            // For m==1, Ez is linear in r, for small r
                            // Therefore, the formula below regularizes the singularity
                            Br(i, j, 0, 2*m-1) += dt*(
                                T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m-1)
                                - m * Ez(i+1, j, 0, 2*m  )/dr );  // Real part
                            Br(i, j, 0, 2*m  ) += dt*(
                                T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m  )
                                + m * Ez(i+1, j, 0, 2*m-1)/dr ); // Imaginary part
                        } else {
                            Br(i, j, 0, 2*m-1) = 0.;
                            Br(i, j, 0, 2*m  ) = 0.;
                        }
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){
                Bt(i, j, 0, 0) += dt*(
                    T_Algo::UpwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 0)
                    - T_Algo::UpwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 0)); // Mode m=0
                for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                    Bt(i, j, 0, 2*m-1) += dt*(
                        T_Algo::UpwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 2*m-1)
                        - T_Algo::UpwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 2*m-1)); // Real part
                    Bt(i, j, 0, 2*m  ) += dt*(
                        T_Algo::UpwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 2*m  )
                        - T_Algo::UpwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 2*m  )); // Imaginary part
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/){
                Real const r = rmin + (i + 0.5)*dr; // r on a cell-centered grid (Bz is cell-centered in r)
                Bz(i, j, 0, 0) += dt*( - T_Algo::UpwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 0));
                for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                    Bz(i, j, 0, 2*m-1) += dt*( m * Er(i, j, 0, 2*m  )/r
                        - T_Algo::UpwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m-1)); // Real part
                    Bz(i, j, 0, 2*m  ) += dt*(-m * Er(i, j, 0, 2*m-1)/r
                        - T_Algo::UpwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m  )); // Imaginary part
                }
            }

        );

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }
}

#endif // corresponds to ifndef WARPX_DIM_RZ
