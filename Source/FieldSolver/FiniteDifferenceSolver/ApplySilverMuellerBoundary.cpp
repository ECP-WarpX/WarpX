/* Copyright 2021
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Utils/WarpXAlgorithmSelection.H"
#include "FiniteDifferenceSolver.H"
#include "Utils/WarpXConst.H"
#include <AMReX_Gpu.H>

using namespace amrex;

/**
 * \brief Update the E field at the boundary, using the Silver-Mueller condition
 */
void FiniteDifferenceSolver::ApplySilverMuellerBoundary (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    amrex::Box domain_box,
    amrex::Real const dt ) {

    // Ensure that we are using the Yee solver
    if (m_fdtd_algo != MaxwellSolverAlgo::Yee) {
        amrex::Abort("The Silver-Mueller boundary conditions can only be used with the Yee solver.");
    }

    // Ensure that we are using the cells the domain
    domain_box.enclosedCells();

    // Calculate relevant coefficients
    amrex::Real const cdt_over_dx = PhysConst::c*dt*m_stencil_coefs_x[0];
    amrex::Real coef1_x = (1. - cdt_over_dx)/(1. + cdt_over_dx);
    amrex::Real coef2_x = 2*cdt_over_dx/(1. + cdt_over_dx) / PhysConst::c;
#ifdef WARPX_DIM_3D
    amrex::Real const cdt_over_dy = PhysConst::c*dt*m_stencil_coefs_y[0];
    amrex::Real coef1_y = (1. - cdt_over_dy)/(1. + cdt_over_dx);
    amrex::Real coef2_y = 2*cdt_over_dy/(1. + cdt_over_dy) / PhysConst::c;
#endif
    amrex::Real const cdt_over_dz = PhysConst::c*dt*m_stencil_coefs_z[0];
    amrex::Real coef1_z = (1. - cdt_over_dz)/(1. + cdt_over_dx);
    amrex::Real coef2_z = 2*cdt_over_dz/(1. + cdt_over_dz) / PhysConst::c;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
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

        // Extract the tileboxes for which to loop
        Box tbx  = mfi.tilebox(Bfield[0]->ixType().toIntVect());
        Box tby  = mfi.tilebox(Bfield[1]->ixType().toIntVect());
        Box tbz  = mfi.tilebox(Bfield[2]->ixType().toIntVect());

        // We will modify the first (i.e. innermost) guard cell
        // (if it is outside of the physical domain)
        // Thus, the tileboxes here are grown by 1 guard cell
        tbx.grow(1);
        tby.grow(1);
        tbz.grow(1);

        // Apply Boundary condition to Bx
        amrex::ParallelFor(tbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k){

#ifdef WARPX_DIM_XZ
                // At the +z boundary (innermost guard cell)
                if ( j==domain_box.bigEnd(1)+1 )
                    Bx(i,j,k) = coef1_z * Bx(i,j,k) - coef2_z * Ey(i,j,k);
                // At the -z boundary (innermost guard cell)
                if ( j==domain_box.smallEnd(1)-1 )
                    Bx(i,j,k) = coef1_z * Bx(i,j,k) + coef2_z * Ey(i,j+1,k);
#endif
            }
        );

        // Apply Boundary condition to By
        amrex::ParallelFor(tby,
            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                // At the +x boundary (innermost guard cell)
                if ( i==domain_box.bigEnd(0)+1 )
                    By(i,j,k) = coef1_x * By(i,j,k) - coef2_x * Ez(i,j,k);
                // At the -x boundary (innermost guard cell)
                if ( i==domain_box.smallEnd(0)-1 )
                    By(i,j,k) = coef1_x * By(i,j,k) + coef2_x * Ez(i+1,j,k);
#ifdef WARPX_DIM_XZ
                // At the +z boundary (innermost guard cell)
                if ( j==domain_box.bigEnd(1)+1 )
                    By(i,j,k) = coef1_z * By(i,j,k) + coef2_z * Ex(i,j,k);
                // At the -z boundary (innermost guard cell)
                if ( j==domain_box.smallEnd(1)-1 )
                    By(i,j,k) = coef1_z * By(i,j,k) - coef2_z * Ex(i,j+1,k);
#endif
            }
        );

        // Apply Boundary condition to Bz
        amrex::ParallelFor(tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k){

                // At the +x boundary (innermost guard cell)
                if ( i==domain_box.bigEnd(0)+1 )
                    Bz(i,j,k) = coef1_x * Bz(i,j,k) + coef2_x * Ey(i,j,k);
                // At the -x boundary (innermost guard cell)
                if ( i==domain_box.smallEnd(0)-1 )
                    Bz(i,j,k) = coef1_x * Bz(i,j,k) - coef2_x * Ey(i+1,j,k);

            }
        );

    }
}
