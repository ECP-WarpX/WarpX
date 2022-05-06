#ifndef NamedComponentParticleContainer_H_
#define NamedComponentParticleContainer_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>

template <int T_NStructReal, int T_NStructInt, int T_NArrayReal, int T_NArrayInt, template<class> class Allocator>
class NamedComponentParticleContainer :
public amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt, Allocator>
{
public:
    NamedComponentParticleContainer () : amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt, Allocator>() {}
    NamedComponentParticleContainer (amrex::AmrParGDB* amr_pgdb) : amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt, Allocator>(amr_pgdb) {}
    virtual ~NamedComponentParticleContainer() {}

    NamedComponentParticleContainer ( NamedComponentParticleContainer && ) = default;
    NamedComponentParticleContainer& operator= ( NamedComponentParticleContainer && ) = default;
};

#endif
