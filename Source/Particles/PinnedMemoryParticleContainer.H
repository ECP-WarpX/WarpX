#ifndef PinnedMemoryParticleContainer_H_
#define PinnedMemoryParticleContainer_H_

#include "NamedComponentParticleContainer.H"
#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <AMReX_ParIter.H>

/*
using PinnedMemoryParticleContainer = amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;
*/

class PinnedMemoryParticleContainer : public NamedComponentParticleContainer<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>
{
public:
    PinnedMemoryParticleContainer () : NamedComponentParticleContainer<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>() {}
    PinnedMemoryParticleContainer (amrex::AmrParGDB* amr_pgdb) : NamedComponentParticleContainer<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>(amr_pgdb) {}
    virtual ~PinnedMemoryParticleContainer() {}

    PinnedMemoryParticleContainer ( PinnedMemoryParticleContainer && ) = default;
    PinnedMemoryParticleContainer& operator= ( PinnedMemoryParticleContainer && ) = default;
};

#endif
