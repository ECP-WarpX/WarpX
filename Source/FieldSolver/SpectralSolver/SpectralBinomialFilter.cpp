/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "SpectralBinomialFilter.H"
#include <AMReX_REAL.H>

#include <cmath>

using namespace amrex::literals;

/* \brief Initialize the general filter array */
void
SpectralBinomialFilter::InitFilterArray (HankelTransform::RealVector const & kvec,
                                         amrex::Real const dels,
                                         int const npasses,
                                         bool const compensation,
                                         KFilterArray & filter)
{

    filter.resize(kvec.size());

    for (int i=0 ; i < kvec.size() ; i++) {
        amrex::Real const ss = std::sin(0.5_rt*kvec[i]*dels);
        amrex::Real const ss2 = ss*ss;
        amrex::Real filt = std::pow(1._rt - ss2, npasses);
        if (compensation) {
            filt *= (1._rt + npasses*ss2);
        }
        filter[i] = filt;
    }

}

/* \brief Initialize the radial filter array */
void
SpectralBinomialFilter::InitFilterArrayR (RealKVector const & kr,
                                          amrex::Real const dr,
                                          int const npasses,
                                          bool const compensation)
{
    // Note that this includes the kr values for all modes
    InitFilterArray(kr, dr, npasses, compensation, filter_r);
}

/* \brief Initialize the longitudinal filter array */
void
SpectralBinomialFilter::InitFilterArrayZ (RealKVector const & kz,
                                          amrex::Real const dz,
                                          int const npasses,
                                          bool const compensation)
{
    InitFilterArray(kz, dz, npasses, compensation, filter_z);
}
