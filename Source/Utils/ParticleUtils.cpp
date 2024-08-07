/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ParticleUtils.H"

#include "WarpX.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Dim3.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

namespace ParticleUtils
{

    using namespace amrex;

    // Define shortcuts for frequently-used type names
    using ParticleType = typename WarpXParticleContainer::ParticleType;
    using ParticleTileType = typename WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = typename ParticleTileType::ParticleTileDataType;
    using ParticleBins = DenseBins<ParticleTileDataType>;
    using index_type = typename ParticleBins::index_type;

    /* Find the particles and count the particles that are in each cell.
       Note that this does *not* rearrange particle arrays */
    ParticleBins
    findParticlesInEachCell (int lev,
                             MFIter const & mfi,
                             ParticleTileType & ptile) {

        // Extract particle structures for this tile
        int const np = ptile.numParticles();
        auto ptd = ptile.getParticleTileData();

        // Extract box properties
        Geometry const& geom = WarpX::GetInstance().Geom(lev);
        Box const& cbx = mfi.tilebox(IntVect::TheZeroVector()); //Cell-centered box
        const auto lo = lbound(cbx);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        // Find particles that are in each cell;
        // results are stored in the object `bins`.
        ParticleBins bins;
        bins.build(np, ptd, cbx,
            // Pass lambda function that returns the cell index
            [=] AMREX_GPU_DEVICE (ParticleType const & p) noexcept -> IntVect
            {
                return getParticleCellIndex(p, plo, dxi, lo);
            });

        return bins;
    }

} // namespace ParticleUtils
