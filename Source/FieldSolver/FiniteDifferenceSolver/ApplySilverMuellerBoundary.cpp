/* Copyright 2021
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

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

    amrex::Print() << "Silver-Mueller: " << domain_box << std::endl;
}
