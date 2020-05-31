/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SpectralKSpaceRZ.H"
#include "SpectralSolverRZ.H"
#include "SpectralAlgorithms/PsatdAlgorithmRZ.H"
#include "SpectralAlgorithms/GalileanPsatdAlgorithmRZ.H"

/* \brief Initialize the spectral Maxwell solver
 *
 * This function selects the spectral algorithm to be used, allocates the
 * corresponding coefficients for the discretized field update equation,
 * and prepares the structures that store the fields in spectral space.
 *
 * \param n_rz_azimuthal_modes Number of azimuthal modes
 * \param norder_z Order of accuracy of the spatial derivatives along z
 * \param nodal    Whether the solver is applied to a nodal or staggered grid
 * \param dx       Cell size along each dimension
 * \param dt       Time step
 * \param pml      Whether the boxes in which the solver is applied are PML boxes
 *                 PML is not supported.
 * \param periodic_single_box Whether the full simulation domain consists of a single periodic box (i.e. the global domain is not MPI parallelized)

 */
SpectralSolverRZ::SpectralSolverRZ(amrex::BoxArray const & realspace_ba,
                                   amrex::DistributionMapping const & dm,
                                   int const n_rz_azimuthal_modes,
                                   int const norder_z, bool const nodal,
                                   const amrex::Array<amrex::Real,3>& v_galilean,
                                   amrex::RealVect const dx, amrex::Real const dt,
                                   int const lev,
                                   bool const pml, bool const periodic_single_box)
{

    // Initialize all structures using the same distribution mapping dm

    // - The k space object contains info about the size of
    //   the spectral space corresponding to each box in `realspace_ba`,
    //   as well as the value of the corresponding k coordinates.
    SpectralKSpaceRZ k_space(realspace_ba, dm, dx);

    // - Select the algorithm depending on the input parameters
    //   Initialize the corresponding coefficients over k space
    //   PML is not supported.
    if (v_galilean[2] == 0) {
         // v_galilean is 0: use standard PSATD algorithm
        algorithm = std::unique_ptr<PsatdAlgorithmRZ>(
            new PsatdAlgorithmRZ(k_space, dm, n_rz_azimuthal_modes, norder_z, nodal, dt));
    } else {
        // Otherwise: use the Galilean algorithm
        algorithm = std::unique_ptr<GalileanPsatdAlgorithmRZ>(
            new GalileanPsatdAlgorithmRZ(k_space, dm, n_rz_azimuthal_modes, norder_z, nodal, v_galilean, dt));
    }

    // - Initialize arrays for fields in spectral space + FFT plans
    field_data = SpectralFieldDataRZ(realspace_ba, k_space, dm,
                                     algorithm->getRequiredNumberOfFields(),
                                     n_rz_azimuthal_modes, lev, periodic_single_box);

};
