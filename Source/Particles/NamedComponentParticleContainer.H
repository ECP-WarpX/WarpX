#ifndef NamedComponentParticleContainer_H_
#define NamedComponentParticleContainer_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>
#include "Utils/TextMsg.H"

struct PIdx
{
    enum { // Particle Attributes stored in amrex::ParticleContainer's struct of array
        w = 0,  // weight
        ux, uy, uz,
#ifdef WARPX_DIM_RZ
theta, // RZ needs all three position components
#endif
nattribs
    };
};

template <template<class> class Allocator>
class NamedComponentParticleContainer :
public amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>
{
public:
    NamedComponentParticleContainer () : amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>() {}
    NamedComponentParticleContainer (amrex::AmrParGDB* amr_pgdb)
    : amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>(amr_pgdb) {

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

    NamedComponentParticleContainer( amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator> && pc )
    : amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>(std::move(pc)) {}

    using amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>::NumRealComps;
    using amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>::NumIntComps;
    using amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>::AddRealComp;
    using amrex::ParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>::AddIntComp;

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

    /** Create an empty particle container
     *
     * This creates a new AMReX particle container type with same compile-time
     * and run-time attributes. But, it can change its allocator. This is
     * helpful when creating temporary particle buffers for filter operations
     * and device-to-host copies.
     */
    template <template<class> class NewAllocator=amrex::DefaultAllocator>
    NamedComponentParticleContainer<NewAllocator>
    make_alike () {
        NamedComponentParticleContainer<NewAllocator> tmp;

        return tmp;
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
