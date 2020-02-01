/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpXAlgorithmSelection.H"
#include "FiniteDifferenceSolver.H"
#ifdef WARPX_DIM_RZ
#   include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
#   include "FiniteDifferenceAlgorithms/YeeAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CKCAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/NodalAlgorithm.H"
#endif
#include <AMReX_Gpu.H>

using namespace amrex;

/**
 * \brief Update the E field, over one timestep
 */
void FiniteDifferenceSolver::EvolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    amrex::Real const dt ) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    if (m_fdtd_algo == MaxwellSolverAlgo::Yee){

        EvolveECylindrical <CylindricalYeeAlgorithm> ( Efield, Bfield, dt );

#else
    if (m_do_nodal) {

        EvolveECartesian <NodalAlgorithm> ( Efield, Bfield, dt );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {

        EvolveECartesian <YeeAlgorithm> ( Efield, Bfield, dt );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        EvolveECartesian <CKCAlgorithm> ( Efield, Bfield, dt );

#endif
    } else {
        amrex::Abort("Unknown algorithm");
    }

}


#ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveECartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    amrex::Real const dt ) {

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        auto const& Ex = Efield[0]->array(mfi);
        auto const& Ey = Efield[1]->array(mfi);
        auto const& Ez = Efield[2]->array(mfi);
        auto const& Bx = Efield[0]->array(mfi);
        auto const& By = Efield[1]->array(mfi);
        auto const& Bz = Efield[2]->array(mfi);

        // Extract stencil coefficients
        Real const* AMREX_RESTRICT coefs_x = stencil_coefs_x.dataPtr();
        int const n_coefs_x = stencil_coefs_x.size();
        Real const* AMREX_RESTRICT coefs_y = stencil_coefs_y.dataPtr();
        int const n_coefs_y = stencil_coefs_y.size();
        Real const* AMREX_RESTRICT coefs_z = stencil_coefs_z.dataPtr();
        int const n_coefs_z = stencil_coefs_z.size();

        // Extract tileboxes for which to loop
        const Box& tex  = mfi.tilebox(Efield[0]->ixType().ixType());
        const Box& tey  = mfi.tilebox(Efield[1]->ixType().ixType());
        const Box& tez  = mfi.tilebox(Efield[2]->ixType().ixType());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tex, tey, tez,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bx(i, j, k) += dt * T_Algo::DownwardDz(Ey, coefs_z, n_coefs_z, i, j, k)
                             - dt * T_Algo::DownwardDy(Ez, coefs_y, n_coefs_y, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                By(i, j, k) += dt * T_Algo::DownwardDx(Ez, coefs_x, n_coefs_x, i, j, k)
                             - dt * T_Algo::DownwardDz(Ex, coefs_z, n_coefs_z, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bz(i, j, k) += dt * T_Algo::DownwardDy(Ex, coefs_y, n_coefs_y, i, j, k)
                             - dt * T_Algo::DownwardDx(Ey, coefs_x, n_coefs_x, i, j, k);
            }

        );

    }

}

#else // corresponds to ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveECylindrical (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    amrex::Real const dt ) {

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        auto const& Er = Efield[0]->array(mfi);
        auto const& Et = Efield[1]->array(mfi);
        auto const& Ez = Efield[2]->array(mfi);
        auto const& Br = Efield[0]->array(mfi);
        auto const& Bt = Efield[1]->array(mfi);
        auto const& Bz = Efield[2]->array(mfi);

        // Extract stencil coefficients
        Real const* AMREX_RESTRICT coefs_r = stencil_coefs_r.dataPtr();
        int const n_coefs_r = stencil_coefs_r.size();
        Real const* AMREX_RESTRICT coefs_z = stencil_coefs_z.dataPtr();
        int const n_coefs_z = stencil_coefs_z.size();

        // Extract cylindrical specific parameters
        Real const dr = m_dr;
        int const nmodes = m_nmodes;
        Real const rmin = m_rmin;

        // Extract tileboxes for which to loop
        const Box& ter  = mfi.tilebox(Efield[0]->ixType().ixType());
        const Box& tet  = mfi.tilebox(Efield[1]->ixType().ixType());
        const Box& tez  = mfi.tilebox(Efield[2]->ixType().ixType());

        // Loop over the cells and update the fields
        amrex::ParallelFor(ter, tet, tez,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Real const r = rmin + i*dr; // r on nodal point (Br is nodal in r)
                Br(i, j, 0, 0) += dt * T_Algo::DownwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 0); // Mode m=0
                for (int m=1; m<nmodes; m++) { // Higher-order modes
                    Br(i, j, 0, 2*m-1) += dt*(
                        T_Algo::DownwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m-1)
                        - m * T_Algo::DivideByR(Ez, r, dr, m, i, j, 0, 2*m  ));  // Real part
                    Br(i, j, 0, 2*m  ) += dt*(
                        T_Algo::DownwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m  )
                        + m * T_Algo::DivideByR(Ez, r, dr, m, i, j, 0, 2*m-1)); // Imaginary part
                }
                // Ensure that Br remains 0 on axis (except for m=1)
                if (r==0) { // On axis
                    Br(i, j, 0, 0) = 0.; // Mode m=0
                    for (int m=2; m<nmodes; m++) { // Higher-order modes (but not m=1)
                        Br(i, j, 0, 2*m-1) = 0.;
                        Br(i, j, 0, 2*m  ) = 0.;
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bt(i, j, 0, 0) += dt*(
                    T_Algo::DownwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 0)
                    - T_Algo::DownwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 0)); // Mode m=0
                for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                    Bt(i, j, 0, 2*m-1) += dt*(
                        T_Algo::DownwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 2*m-1)
                        - T_Algo::DownwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 2*m-1)); // Real part
                    Bt(i, j, 0, 2*m  ) += dt*(
                        T_Algo::DownwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 2*m  )
                        - T_Algo::DownwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 2*m  )); // Imaginary part
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Real const r = rmin + (i + 0.5)*dr; // r on a cell-centered grid (Bz is cell-centered in r)
                Bz(i, j, 0, 0) += dt*( - T_Algo::DownwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 0));
                for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                    Bz(i, j, 0, 2*m-1) += dt*( m * Er(i, j, 0, 2*m  )/r
                        - T_Algo::DownwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m-1)); // Real part
                    Bz(i, j, 0, 2*m  ) += dt*(-m * Er(i, j, 0, 2*m-1)/r
                        - T_Algo::DownwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m  )); // Imaginary part
                }
            }

        );

    }

}

#endif // corresponds to ifndef WARPX_DIM_RZ
