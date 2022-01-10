#ifndef PinnedMemoryParticleContainer_H_
#define PinnedMemoryParticleContainer_H_

#include "WarpXParticleContainer.H"
#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <AMReX_ParIter.H>

using PinnedMemoryParticleContainer = typename amrex::AmrParticleContainer<0,0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;

#endif
