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
#include "WarpX.H"
#include "Utils/WarpXProfilerWrapper.H"

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
 */
SpectralSolverRZ::SpectralSolverRZ (amrex::BoxArray const & realspace_ba,
                                    amrex::DistributionMapping const & dm,
                                    int const n_rz_azimuthal_modes,
                                    int const norder_z, bool const nodal,
                                    const amrex::Array<amrex::Real,3>& v_galilean,
                                    amrex::RealVect const dx, amrex::Real const dt,
                                    int const lev)
    : k_space(realspace_ba, dm, dx)
{
    // Initialize all structures using the same distribution mapping dm

    // - The k space object contains info about the size of
    //   the spectral space corresponding to each box in `realspace_ba`,
    //   as well as the value of the corresponding k coordinates.

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
                                     n_rz_azimuthal_modes, lev);
}

/* \brief Transform the component `i_comp` of MultiFab `field_mf`
 *  to spectral space, and store the corresponding result internally
 *  (in the spectral field specified by `field_index`) */
void
SpectralSolverRZ::ForwardTransform (amrex::MultiFab const & field_mf, int const field_index,
                                    int const i_comp) {
    WARPX_PROFILE("SpectralSolverRZ::ForwardTransform");
    field_data.ForwardTransform(field_mf, field_index, i_comp);
}

/* \brief Transform the two MultiFabs `field_mf1` and `field_mf2`
 *  to spectral space, and store the corresponding results internally
 *  (in the spectral field specified by `field_index1` and `field_index2`) */
void
SpectralSolverRZ::ForwardTransform (amrex::MultiFab const & field_mf1, int const field_index1,
                                    amrex::MultiFab const & field_mf2, int const field_index2) {
    WARPX_PROFILE("SpectralSolverRZ::ForwardTransform");
    field_data.ForwardTransform(field_mf1, field_index1,
                                field_mf2, field_index2);
}

/* \brief Transform spectral field specified by `field_index` back to
 * real space, and store it in the component `i_comp` of `field_mf` */
void
SpectralSolverRZ::BackwardTransform (amrex::MultiFab& field_mf, int const field_index,
                                     int const i_comp) {
    WARPX_PROFILE("SpectralSolverRZ::BackwardTransform");
    field_data.BackwardTransform(field_mf, field_index, i_comp);
}

/* \brief Transform spectral fields specified by `field_index1` and `field_index2`
 * back to real space, and store it in `field_mf1` and `field_mf2`*/
void
SpectralSolverRZ::BackwardTransform (amrex::MultiFab& field_mf1, int const field_index1,
                                     amrex::MultiFab& field_mf2, int const field_index2) {
    WARPX_PROFILE("SpectralSolverRZ::BackwardTransform");
    field_data.BackwardTransform(field_mf1, field_index1,
                                 field_mf2, field_index2);
}

/* \brief Update the fields in spectral space, over one timestep */
void
SpectralSolverRZ::pushSpectralFields () {
    WARPX_PROFILE("SpectralSolverRZ::pushSpectralFields");
    // Virtual function: the actual function used here depends
    // on the sub-class of `SpectralBaseAlgorithm` that was
    // initialized in the constructor of `SpectralSolverRZ`
    algorithm->pushSpectralFields(field_data);
}

/**
  * \brief Public interface to call the member function ComputeSpectralDivE
  * of the base class SpectralBaseAlgorithmRZ from objects of class SpectralSolverRZ
  */
void
SpectralSolverRZ::ComputeSpectralDivE (const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                                       amrex::MultiFab& divE) {
    algorithm->ComputeSpectralDivE(field_data, Efield, divE);
}

/**
 * \brief Public interface to call the virtual function \c CurrentCorrection,
 * defined in the base class SpectralBaseAlgorithmRZ and possibly overridden
 * by its derived classes (e.g. PsatdAlgorithmRZ), from
 * objects of class SpectralSolverRZ through the private unique pointer \c algorithm
 *
 * \param[in,out] current two-dimensional array of unique pointers to MultiFab
 *                        storing the three components of the current density
 * \param[in]     rho     unique pointer to MultiFab storing the charge density
 */
void
SpectralSolverRZ::CurrentCorrection (std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                      const std::unique_ptr<amrex::MultiFab>& rho) {
     algorithm->CurrentCorrection(field_data, current, rho);
}

void
SpectralSolverRZ::VayDeposition (std::array<std::unique_ptr<amrex::MultiFab>,3>& current)
{
  algorithm->VayDeposition(field_data, current);
}
