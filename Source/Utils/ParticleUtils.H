/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLE_UTILS_H_
#define WARPX_PARTICLE_UTILS_H_

#include "Particles/WarpXParticleContainer.H"
#include <AMReX_DenseBins.H>

namespace ParticleUtils {

    /**
    * \brief Find the particles and count the particles that are in each cell. More specifically
    * this function returns an amrex::DenseBins object containing an offset array and a permutation
    * array which can be used to loop over all the cells in a tile and apply an algorithm to
    * particles of a given species present in each cell.
    * Note that this does *not* rearrange particle arrays.
    *
    * @param[in] lev the index of the refinement level.
    * @param[in] mfi the MultiFAB iterator.
    * @param[in] ptile the particle tile.
    */
    amrex::DenseBins<WarpXParticleContainer::ParticleType>
    findParticlesInEachCell( int const lev, amrex::MFIter const& mfi,
                             WarpXParticleContainer::ParticleTileType const& ptile);
}

#endif // WARPX_PARTICLE_UTILS_H_
