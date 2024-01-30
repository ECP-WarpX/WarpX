/* Copyright 2021-2023 The WarpX Community
 *
 * Authors: Axel Huebl, Remi Lehe, Roelof Groenewald
 * License: BSD-3-Clause-LBNL
 */

#include "Python/pyWarpX.H"

#include <Particles/PinnedMemoryParticleContainer.H>


void init_PinnedMemoryParticleContainer (py::module& m)
{
    py::class_<
        PinnedMemoryParticleContainer,
        amrex::ParticleContainer<0,0,PIdx::nattribs,0,amrex::PinnedArenaAllocator>
    > pmpc (m, "PinnedMemoryParticleContainer");
}
