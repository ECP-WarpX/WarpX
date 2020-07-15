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

/* \brief Initialize the radial filter array */
void
SpectralBinomialFilter::InitFilterArrayR (HankelTransform::RealVector const & kr,
                                          amrex::Real const dr,
                                          int const npasses,
                                          bool const compensation)
{

    // Note that this includes the kr values for all modes
    filter_r.resize(kr.size());

    for (int i=0 ; i < kr.size() ; i++) {
        amrex::Real const sr = std::sin(0.5_rt*kr[i]*dr);
        amrex::Real const sr2 = sr*sr;
        amrex::Real filt_r = std::pow(1._rt - sr2, npasses);
        if (compensation) {
            filt_r *= (1._rt + npasses*sr2);
        }
        filter_r[i] = filt_r;
    }

}

/* \brief Initialize the longitudinal filter array */
void
SpectralBinomialFilter::InitFilterArrayZ (RealKVector const & kz,
                                          amrex::Real const dz,
                                          int const npasses,
                                          bool const compensation)
{

    filter_z.resize(kz.size());

    for (int i=0 ; i < kz.size() ; i++) {
        amrex::Real const sz = std::sin(0.5_rt*kz[i]*dz);
        amrex::Real const sz2 = sz*sz;
        amrex::Real filt_z = std::pow(1._rt - sz2, npasses);
        if (compensation) {
            filt_z *= (1._rt + npasses*sz2);
        }
        filter_z[i] = filt_z;
    }
}
