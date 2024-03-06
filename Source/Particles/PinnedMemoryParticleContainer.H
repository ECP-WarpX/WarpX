#ifndef PinnedMemoryParticleContainer_H_
#define PinnedMemoryParticleContainer_H_

#include "NamedComponentParticleContainer.H"

using PinnedMemoryParticleContainer = NamedComponentParticleContainer<amrex::PinnedArenaAllocator>;

#endif
