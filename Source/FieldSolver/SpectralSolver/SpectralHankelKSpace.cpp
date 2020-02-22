/* Copyright 2019 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <WarpXConst.H>
#include <SpectralHankelKSpace.H>
#include <cmath>

using namespace amrex;

/* \brief Initialize k space object.
 *
 * \param realspace_ba Box array that corresponds to the decomposition
 * of the fields in real space (cell-centered ; includes guard cells)
 * \param dm Indicates which MPI proc owns which box, in realspace_ba.
 * \param realspace_dx Cell size of the grid in real space
 */
SpectralHankelKSpace::SpectralHankelKSpace (const BoxArray& realspace_ba,
                                            const DistributionMapping& dm,
                                            const RealVect realspace_dx)
{
    dx = realspace_dx;  // Store the cell size as member `dx`

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        realspace_ba.ixType() == IndexType::TheCellType(),
        "SpectralHankelKSpace expects a cell-centered box.");

    // Create the box array that corresponds to spectral space
    BoxList spectral_bl; // Create empty box list

    // Loop over boxes and fill the box list
    for (int i=0; i < realspace_ba.size(); i++ ) {
        // For local FFTs, boxes in spectral space start at 0 in
        // each direction and have the same number of points as the
        // (cell-centered) real space box
        Box realspace_bx = realspace_ba[i];
        IntVect fft_size = realspace_bx.length();
        IntVect spectral_bx_size = fft_size;
        // Define the corresponding box
        Box spectral_bx = Box(IntVect::TheZeroVector(),
                              spectral_bx_size - IntVect::TheUnitVector() );
        spectral_bl.push_back(spectral_bx);
    }
    spectralspace_ba.define(spectral_bl);

    // Allocate the components of the kz vector
    const int i_dim = 1;
    const bool only_positive_k = false;
    k_vec[i_dim] = getKComponent(dm, realspace_ba, i_dim, only_positive_k);

}
