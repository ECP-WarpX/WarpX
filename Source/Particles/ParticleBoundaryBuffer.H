/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARYBUFFER_H_
#define PARTICLEBOUNDARYBUFFER_H_

#include "Particles/MultiParticleContainer_fwd.H"
#include "Particles/WarpXParticleContainer.H"
#include "Particles/PinnedMemoryParticleContainer.H"

#include <vector>

/**
 *  This stores particles that have left / been absorbed by domain and embedded boundaries.
 */
class ParticleBoundaryBuffer
{
public:
    ParticleBoundaryBuffer ();

    int numSpecies() const { return getSpeciesNames().size(); }

    const std::vector<std::string>& getSpeciesNames() const {
        static bool initialized = false;
        if (!initialized)
        {
            amrex::ParmParse pp_particles("particles");
            pp_particles.queryarr("species_names", m_species_names);
            initialized = true;
        }
        return m_species_names;
    }

    void gatherParticles (MultiParticleContainer& mypc,
                          const amrex::Vector<const amrex::MultiFab*>& distance_to_eb);

    void redistribute ();
    void clearParticles ();
    void clearParticles (int const i);

    void printNumParticles () const;

    int getNumParticlesInContainer(const std::string species_name, int boundary);

    PinnedMemoryParticleContainer& getParticleBuffer(const std::string species_name, int boundary);

    PinnedMemoryParticleContainer* getParticleBufferPointer(const std::string species_name, int boundary);

    static constexpr int numBoundaries () {
        return AMREX_SPACEDIM*2
#ifdef AMREX_USE_EB
            + 1
#endif
            ;
    }

    bool isDefinedForAnySpecies (int const ibuffer) {return (m_do_any_boundary[ibuffer] != 0);}

    std::string boundaryName (int const ibuffer) {return m_boundary_names[ibuffer];}

private:
    // over boundary, then number of species
    std::vector<std::vector<PinnedMemoryParticleContainer> > m_particle_containers;

    // over boundary, then number of species
    std::vector<std::vector<int> > m_do_boundary_buffer;

    // over boundary
    std::vector<int> m_do_any_boundary;
    std::vector<std::string> m_boundary_names;

    mutable std::vector<std::string> m_species_names;
};

#endif /*PARTICLEBOUNDARYBUFFER_H_*/
