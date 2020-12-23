/* Copyright 2019 Axel Huebl, David Grote, Igor Andriyash
 * Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FielIO_H_
#define WARPX_FielIO_H_

#include <WarpX.H>

void
AverageAndPackVectorField( amrex::MultiFab& mf_avg,
                         const std::array< std::unique_ptr<amrex::MultiFab>, 3 >& vector_field,
                         const amrex::DistributionMapping& dm,
                         const int dcomp, const int ngrow );

void
AverageAndPackScalarField( amrex::MultiFab& mf_avg,
                         const amrex::MultiFab & scalar_field,
                         const amrex::DistributionMapping& dm,
                         const int dcomp, const int ngrow );

std::vector<std::uint64_t>
getReversedVec( const amrex::IntVect& v );

std::vector<double>
getReversedVec( const amrex::Real* v );

#endif // WARPX_FielIO_H_
