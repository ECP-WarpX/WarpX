/* Copyright 2023 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_GENERICPARTICLEFUNCTOR_H_
#define WARPX_GENERICPARTICLEFUNCTOR_H_

#include "ComputeParticleDiagFunctor.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <AMReX.H>
#include <AMReX_AmrParticles.H>

/**
 * \brief Functor to select particles and store them in particle buffers
 */
template <typename PC>
class
GenericParticleFunctor final : public ComputeParticleDiagFunctor
{
public:
    /** Selects particles to be written and copies them to pc_dst*/
    GenericParticleFunctor(PC* pc_src, std::string species_name): m_pc_src{pc_src}, m_species_name{std::move(species_name)} {}
    void operator () (PinnedMemoryParticleContainer& pc_dst, int &TotalParticleCounter, int i_buffer) const override {}

private:
    /** Source particle container */
    PC* m_pc_src = nullptr;
    /** String containing species name of particles being written */
    std::string m_species_name;
};

#endif //
