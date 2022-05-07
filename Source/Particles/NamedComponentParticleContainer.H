#ifndef NamedComponentParticleContainer_H_
#define NamedComponentParticleContainer_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>
#include "Utils/TextMsg.H"

template <int T_NStructReal, int T_NStructInt, int T_NArrayReal, int T_NArrayInt, template<class> class Allocator>
class NamedComponentParticleContainer :
public amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt, Allocator>
{
public:
    NamedComponentParticleContainer () : amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt, Allocator>() {}
    NamedComponentParticleContainer (amrex::AmrParGDB* amr_pgdb)
    : amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt, Allocator>(amr_pgdb) {

        // build up the map of string names to particle component numbers
        particle_comps["w"]  = PIdx::w;
        particle_comps["ux"] = PIdx::ux;
        particle_comps["uy"] = PIdx::uy;
        particle_comps["uz"] = PIdx::uz;
    #ifdef WARPX_DIM_RZ
        particle_comps["theta"] = PIdx::theta;
    #endif

    }
    virtual ~NamedComponentParticleContainer() {}

    NamedComponentParticleContainer ( NamedComponentParticleContainer && ) = default;
    NamedComponentParticleContainer& operator= ( NamedComponentParticleContainer && ) = default;

    using amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt>::NumRealComps;
    using amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt>::NumIntComps;
    using amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt>::AddRealComp;
    using amrex::ParticleContainer<T_NStructReal, T_NStructInt, T_NArrayReal, T_NArrayInt>::AddIntComp;

    void AddRealComp (const std::string& name, bool comm=true)
    {
        auto search = particle_comps.find(name);
        if (search == particle_comps.end()) {
            particle_comps[name] = NumRealComps();
            particle_runtime_comps[name] = NumRealComps() - PIdx::nattribs;
            AddRealComp(comm);
        } else {
            amrex::Print() << Utils::TextMsg::Info(
                name + " already exists in particle_comps, not adding.");
        }
    }

    void AddIntComp (const std::string& name, bool comm=true)
    {
        auto search = particle_icomps.find(name);
        if (search == particle_icomps.end()) {
            particle_icomps[name] = NumIntComps();
            particle_runtime_icomps[name] = NumIntComps() - 0;
            AddIntComp(comm);
        } else {
            amrex::Print() << Utils::TextMsg::Info(
                name + " already exists in particle_icomps, not adding.");
        }
    }

    std::map<std::string, int> getParticleComps () const noexcept { return particle_comps;}
    std::map<std::string, int> getParticleiComps () const noexcept { return particle_icomps;}
    std::map<std::string, int> getParticleRuntimeComps () const noexcept { return particle_runtime_comps;}
    std::map<std::string, int> getParticleRuntimeiComps () const noexcept { return particle_runtime_icomps;}

protected:
    std::map<std::string, int> particle_comps;
    std::map<std::string, int> particle_icomps;
    std::map<std::string, int> particle_runtime_comps;
    std::map<std::string, int> particle_runtime_icomps;
};

#endif
