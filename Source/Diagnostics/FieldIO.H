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
#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif

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

#ifdef WARPX_USE_OPENPMD
void
setOpenPMDUnit( openPMD::Mesh mesh, const std::string field_name );

std::vector<std::uint64_t>
getReversedVec( const amrex::IntVect& v );

std::vector<double>
getReversedVec( const amrex::Real* v );

#endif // WARPX_USE_OPENPMD

#endif // WARPX_FielIO_H_
