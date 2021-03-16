/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SpectralKSpace.H"
#include "SpectralSolver.H"
#include "SpectralAlgorithms/PsatdAlgorithm.H"
#include "SpectralAlgorithms/PMLPsatdAlgorithm.H"
#include "SpectralAlgorithms/ComovingPsatdAlgorithm.H"
#include "WarpX.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"

#include <memory>

#if WARPX_USE_PSATD

/* \brief Initialize the spectral Maxwell solver
 *
 * This function selects the spectral algorithm to be used, allocates the
 * corresponding coefficients for the discretized field update equation,
 * and prepares the structures that store the fields in spectral space.
 *
 * \param norder_x Order of accuracy of the spatial derivatives along x
 * \param norder_y Order of accuracy of the spatial derivatives along y
 * \param norder_z Order of accuracy of the spatial derivatives along z
 * \param nodal    Whether the solver is applied to a nodal or staggered grid
 * \param dx       Cell size along each dimension
 * \param dt       Time step
 * \param pml      Whether the boxes in which the solver is applied are PML boxes
 * \param periodic_single_box Whether the full simulation domain consists of a single periodic box (i.e. the global domain is not MPI parallelized)
 */
SpectralSolver::SpectralSolver(
                const int lev,
                const amrex::BoxArray& realspace_ba,
                const amrex::DistributionMapping& dm,
                const int norder_x, const int norder_y,
                const int norder_z, const bool nodal,
                const amrex::Array<amrex::Real,3>& v_galilean,
                const amrex::Array<amrex::Real,3>& v_comoving,
                const amrex::RealVect dx, const amrex::Real dt,
                const bool pml, const bool periodic_single_box,
                const bool update_with_rho,
                const bool fft_do_time_averaging) {

    // Initialize all structures using the same distribution mapping dm

    // - Initialize k space object (Contains info about the size of
    // the spectral space corresponding to each box in `realspace_ba`,
    // as well as the value of the corresponding k coordinates)
    const SpectralKSpace k_space= SpectralKSpace(realspace_ba, dm, dx);

    // - Select the algorithm depending on the input parameters
    //   Initialize the corresponding coefficients over k space

    if (pml) {
        algorithm = std::make_unique<PMLPsatdAlgorithm>(
            k_space, dm, norder_x, norder_y, norder_z, nodal, dt);
    }
    else {
        // Comoving PSATD algorithm
        if (v_comoving[0] != 0. || v_comoving[1] != 0. || v_comoving[2] != 0.) {
            algorithm = std::make_unique<ComovingPsatdAlgorithm>(
                k_space, dm, norder_x, norder_y, norder_z, nodal, v_comoving, dt, update_with_rho);
        }
        // PSATD algorithms: standard, Galilean, or averaged Galilean
        else {
            algorithm = std::make_unique<PsatdAlgorithm>(
                k_space, dm, norder_x, norder_y, norder_z, nodal, v_galilean, dt, update_with_rho, fft_do_time_averaging);
        }
    }

    // - Initialize arrays for fields in spectral space + FFT plans
    field_data = SpectralFieldData( lev, realspace_ba, k_space, dm,
                    algorithm->getRequiredNumberOfFields(), periodic_single_box);

}

void
SpectralSolver::ForwardTransform( const int lev,
                                  const amrex::MultiFab& mf,
                                  const int field_index,
                                  const int i_comp )
{
    WARPX_PROFILE("SpectralSolver::ForwardTransform");
    field_data.ForwardTransform( lev, mf, field_index, i_comp );
}

void
SpectralSolver::BackwardTransform( const int lev,
                                   amrex::MultiFab& mf,
                                   const int field_index,
                                   const int i_comp )
{
    WARPX_PROFILE("SpectralSolver::BackwardTransform");
    field_data.BackwardTransform( lev, mf, field_index, i_comp );
}

void
SpectralSolver::pushSpectralFields(){
    WARPX_PROFILE("SpectralSolver::pushSpectralFields");
    // Virtual function: the actual function used here depends
    // on the sub-class of `SpectralBaseAlgorithm` that was
    // initialized in the constructor of `SpectralSolver`
    algorithm->pushSpectralFields( field_data );
}

#endif // WARPX_USE_PSATD
