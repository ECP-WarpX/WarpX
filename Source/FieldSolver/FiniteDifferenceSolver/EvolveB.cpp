/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpX.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "FiniteDifferenceSolver.H"
#ifdef WARPX_DIM_RZ
#   include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#endif
#include <AMReX_Gpu.H>

using namespace amrex;

/**
 * \brief Update the B field, over one timestep
 */
void FiniteDifferenceSolver::EvolveB (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& face_areas,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& area_enl,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& area_red,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& area_stab,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Vfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Rhofield,
    std::array< std::unique_ptr<amrex::iMultiFab>, 3 >& flag_unst_cell,
    std::array< std::unique_ptr<amrex::LayoutData<FaceInfoBox> >, 3 >& borrowing,
    std::array< std::unique_ptr<amrex::LayoutData<FaceInfoBox> >, 3 >& lending,
    int lev, amrex::Real const dt ) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    if (m_fdtd_algo == MaxwellSolverAlgo::Yee){

        EvolveBCylindrical <CylindricalYeeAlgorithm> ( Bfield, Efield, lev, dt );

#else
    if (m_do_nodal) {

        EvolveBCartesian <CartesianNodalAlgorithm> ( Bfield, Efield, face_areas, lev, dt );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {

        EvolveBCartesian <CartesianYeeAlgorithm> ( Bfield, Efield, face_areas, lev, dt );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        EvolveBCartesian <CartesianCKCAlgorithm> ( Bfield, Efield, face_areas, lev, dt );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::ECT) {

        EvolveRhoCartesianECT( Efield, edge_lengths, face_areas, Rhofield, lev );
        EvolveBCartesianECT( Bfield, face_areas, area_enl, area_red, area_stab, Rhofield,
                             flag_unst_cell, borrowing, lending, lev, dt);
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
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& face_areas,
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

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }
}

void FiniteDifferenceSolver::EvolveRhoCartesianECT (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& face_areas,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Rhofield, int lev ) {

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Rhofield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        // Extract field data for this grid/tile
        Array4<Real> const &Ex = Efield[0]->array(mfi);
        Array4<Real> const &Ey = Efield[1]->array(mfi);
        Array4<Real> const &Ez = Efield[2]->array(mfi);
        Array4<Real> const &Rhox = Rhofield[0]->array(mfi);
        Array4<Real> const &Rhoy = Rhofield[1]->array(mfi);
        Array4<Real> const &Rhoz = Rhofield[2]->array(mfi);
        amrex::Array4<amrex::Real> const &lx = edge_lengths[0]->array(mfi);
        amrex::Array4<amrex::Real> const &ly = edge_lengths[1]->array(mfi);
        amrex::Array4<amrex::Real> const &lz = edge_lengths[2]->array(mfi);
        amrex::Array4<amrex::Real> const &Sx = face_areas[0]->array(mfi);
        amrex::Array4<amrex::Real> const &Sy = face_areas[1]->array(mfi);
        amrex::Array4<amrex::Real> const &Sz = face_areas[2]->array(mfi);

        // Extract tileboxes for which to loop
        Box const &trhox = mfi.tilebox(Rhofield[0]->ixType().toIntVect());
        Box const &trhoy = mfi.tilebox(Rhofield[1]->ixType().toIntVect());
        Box const &trhoz = mfi.tilebox(Rhofield[2]->ixType().toIntVect());

        amrex::ParallelFor(trhox, trhoy, trhoz,

            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (Sx(i, j, k) <= 0 or isnan(Sx(i, j, k))) return;

                Rhox(i, j, k) = (Ey(i, j, k) * ly(i, j, k) - Ey(i, j, k + 1) * ly(i, j, k + 1) +
                    Ez(i, j + 1, k) * lz(i, j + 1, k) - Ez(i, j, k) * lz(i, j, k)) / Sx(i, j, k);

            },

            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (Sy(i, j, k) <= 0 or isnan(Sy(i, j, k))) return;

                Rhoy(i, j, k) = (Ez(i, j, k) * lz(i, j, k) - Ez(i + 1, j, k) * lz(i + 1, j, k) +
                    Ex(i, j, k + 1) * lx(i, j, k + 1) - Ex(i, j, k) * lx(i, j, k)) / Sy(i, j, k);

            },

            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (Sz(i, j, k) <= 0 or isnan(Sz(i, j, k))) return;

                Rhoz(i, j, k) =  (Ex(i, j, k) * lx(i, j, k) - Ex(i, j + 1, k) * lx(i, j + 1, k) +
                    Ey(i + 1, j, k) * ly(i + 1, j, k) - Ey(i, j, k) * ly(i, j, k)) / Sz(i, j, k);

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

void FiniteDifferenceSolver::EvolveBCartesianECT (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& face_areas,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& area_enl,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& area_red,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& area_stab,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Rhofield,
    std::array< std::unique_ptr<amrex::iMultiFab>, 3 >& flag_unst_cell,
    std::array< std::unique_ptr<amrex::LayoutData<FaceInfoBox> >, 3 >& borrowing,
    std::array< std::unique_ptr<amrex::LayoutData<FaceInfoBox> >, 3 >& lending,
    int lev, amrex::Real const dt ) {

    amrex::LayoutData<amrex::Real> *cost = WarpX::getCosts(lev);

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*Bfield[0]); mfi.isValid(); ++mfi) {

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Extract field data for this grid/tile
            Array4<Real> const &B = Bfield[idim]->array(mfi);
            Array4<Real> const &Rho = Rhofield[idim]->array(mfi);
            amrex::Array4<int> const &flag_unst_cell_dim = flag_unst_cell[idim]->array(mfi);
            amrex::Array4<Real> const &S = face_areas[idim]->array(mfi);
            amrex::Array4<Real> const &S_enl = area_enl[idim]->array(mfi);
            amrex::Array4<Real> const &S_red = area_red[idim]->array(mfi);
            amrex::Array4<Real> const &S_stab = area_stab[idim]->array(mfi);

            auto &lending_dim = (*lending[idim])[mfi];
            auto &borrowing_dim = (*borrowing[idim])[mfi];

            auto const &borrowing_inds = (*borrowing[idim])[mfi].inds.array();
            auto const &lending_inds = (*lending[idim])[mfi].inds.array();

            // Extract tileboxes for which to loop
            Box const &tb = mfi.tilebox(Bfield[idim]->ixType().toIntVect());

            //Take care of the unstable cells
            amrex::LoopOnCpu(tb,

            [=, &borrowing_dim, &lending_dim] (int i, int j, int k) {

                if (S(i, j, k) <= 0 or isnan(S(i, j, k))) return;

                if (!flag_unst_cell_dim(i, j, k))
                    return;

                amrex::Real V_enl = Rho(i, j, k) * S(i, j, k);
                amrex::Real rho_enl;
                if (borrowing_inds(i, j, k).size() == 0) {
                    amrex::Abort("EvolveBCartesianECT: face ("
                                + std::to_string(i) + ", "
                                + std::to_string(j) + ", "
                                + std::to_string(k)
                                + ") wasn't extended correctly");
                }
                // First we compute the rho of the enlarged face
                for (int ind : borrowing_inds(i, j, k)) {
                    int ip = borrowing_dim.i_face[ind];
                    int jp = borrowing_dim.j_face[ind];
                    int kp = borrowing_dim.k_face[ind];
                    V_enl += Rho(ip, jp, kp) * borrowing_dim.area[ind];
                }

                rho_enl = V_enl / S_enl(i, j, k);

                //Now we have to insert the computed rho_enl in the lending FaceInfoBox in the correct
                // position
                for (int ind : borrowing_inds(i, j, k)) {
                    int ip = borrowing_dim.i_face[ind];
                    int jp = borrowing_dim.j_face[ind];
                    int kp = borrowing_dim.k_face[ind];
                    for (int ind2 : lending_inds(ip, jp, kp)) {
                        int ip2 = lending_dim.i_face[ind2];
                        int jp2 = lending_dim.j_face[ind2];
                        int kp2 = lending_dim.k_face[ind2];
                        if (ip2 == i and jp2 == j and kp2 == k) {
                            lending_dim.rho_face[ind2] = rho_enl;
                        }
                    }
                }

                B(i, j, k) = B(i, j, k) - dt * rho_enl;

            });

            //Take care of the stable cells
            amrex::ParallelFor(tb,

            [=, &lending_dim] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (S(i, j, k) <= 0 or isnan(S(i, j, k))) return;

                if (flag_unst_cell_dim(i, j, k))
                    return;

                if (lending_inds(i, j, k).size()  == 0) {
                    //Stable cell which hasn't been intruded
                    B(i, j, k) = B(i, j, k) - dt * Rho(i, j, k);
                } else {

                    //Stable cell which has been intruded
                    amrex::Real Venl = Rho(i, j, k) * S_red(i, j, k);

                    for (int ind : lending_inds(i, j, k)) {
                        Venl += lending_dim.rho_face[ind] * lending_dim.area[ind];
                    }

                    B(i, j, k) = B(i, j, k) - dt * Venl / S(i, j, k);
                }

            });

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
