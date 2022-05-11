#ifndef NamedComponentParticleContainer_H_
#define NamedComponentParticleContainer_H_

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
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

/** Particle Container class that allows to add/access particle components
 *  with a name (string) instead of doing so with an integer index.
 *
 *  This is done by storing maps that give the index of the component
 *  that corresponds to a given string.
 */
template <template<class> class Allocator>
class NamedComponentParticleContainer :
public amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>
{
public:
    NamedComponentParticleContainer () : amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>() {}
    NamedComponentParticleContainer (amrex::AmrParGDB* amr_pgdb)
    : amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>(amr_pgdb) {
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

    /** Construct a NamedComponentParticleContainer from a regular
      * amrex::ParticleContainer, and additional name-to-index maps */
    NamedComponentParticleContainer(
        amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator> && pc,
        std::map<std::string, int> p_comps,
        std::map<std::string, int> p_icomps,
        std::map<std::string, int> p_rcomps,
        std::map<std::string, int> p_ricomps)
    : amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>(std::move(pc)),
    particle_comps(p_comps),
    particle_icomps(p_icomps),
    particle_runtime_comps(p_rcomps),
    particle_runtime_icomps(p_ricomps) {}

    /* Move operators */
    NamedComponentParticleContainer ( NamedComponentParticleContainer && ) = default;
    NamedComponentParticleContainer& operator= ( NamedComponentParticleContainer && ) = default;

    /** Create an empty particle container
     *
     * This creates a new NamedComponentParticleContainer with same compile-time
     * and run-time attributes. But it can change its allocator.
     *
     * This function overloads the corresponding function from the parent
     * class (amrex::ParticleContainer)
     */
    template <template<class> class NewAllocator=amrex::DefaultAllocator>
    NamedComponentParticleContainer<NewAllocator>
    make_alike () const {
        auto tmp = NamedComponentParticleContainer<NewAllocator>(
            amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>::template make_alike<NewAllocator>(),
            particle_comps,
            particle_icomps,
            particle_runtime_comps,
            particle_runtime_icomps);

        return tmp;
    }

    using amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>::NumRealComps;
    using amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>::NumIntComps;
    using amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>::AddRealComp;
    using amrex::ParticleContainer<0,0,PIdx::nattribs,0,Allocator>::AddIntComp;

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
