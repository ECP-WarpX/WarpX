/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ParticleUtils.H"
#include "WarpX.H"

namespace ParticleUtils {

    using namespace amrex;
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleBins = DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;

    /* Find the particles and count the particles that are in each cell.
       Note that this does *not* rearrange particle arrays */
    ParticleBins
    findParticlesInEachCell( int const lev, MFIter const& mfi,
                             ParticleTileType const& ptile) {

        // Extract particle structures for this tile
        int const np = ptile.numParticles();
        ParticleType const* particle_ptr = ptile.GetArrayOfStructs()().data();

        // Extract box properties
        Geometry const& geom = WarpX::GetInstance().Geom(lev);
        Box const& cbx = mfi.tilebox(IntVect::TheZeroVector()); //Cell-centered box
        const auto lo = lbound(cbx);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        // Find particles that are in each cell;
        // results are stored in the object `bins`.
        ParticleBins bins;
        bins.build(np, particle_ptr, cbx,
            // Pass lambda function that returns the cell index
            [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> IntVect
            {
                return IntVect(AMREX_D_DECL(
                                   static_cast<int>((p.pos(0)-plo[0])*dxi[0] - lo.x),
                                   static_cast<int>((p.pos(1)-plo[1])*dxi[1] - lo.y),
                                   static_cast<int>((p.pos(2)-plo[2])*dxi[2] - lo.z)));
            });

        return bins;
    }

}
