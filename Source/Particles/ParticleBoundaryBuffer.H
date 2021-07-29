/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARYBUFFER_H_
#define PARTICLEBOUNDARYBUFFER_H_

#include "Particles/ParticleBuffer.H"
#include "Particles/MultiParticleContainer_fwd.H"

#include <vector>

/**
 *  This stores particles that have left / been absorbed by domain and embedded boundaries.
 */
class ParticleBoundaryBuffer
{
public:
    ParticleBoundaryBuffer ();

    int nSpecies() const {return m_species_names.size();}
    std::vector<std::string> GetSpeciesNames() const { return m_species_names; }

    void gatherParticles (MultiParticleContainer& mypc);

    void printNumParticles () const;

private:
    // over boundary, then number of species
    std::vector<std::vector<ParticleBuffer::BufferType<amrex::PinnedArenaAllocator> > > m_particle_containers;

    // over boundary, then number of species
    std::vector<std::vector<int> > m_do_boundary_buffer;

    std::vector<std::string> m_species_names;
};

#endif /*PARTICLEBOUNDARYBUFFER_H_*/
