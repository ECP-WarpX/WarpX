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
class
GenericParticleFunctor final : public ComputeParticleDiagFunctor
{
public:
    GenericParticleFunctor(WarpXParticleContainer *pc_src, std::string species_name);
    /** Selects particles to be written and copies them to pc_dst*/
    void operator () (PinnedMemoryParticleContainer& pc_dst, int &TotalParticleCounter, int i_buffer) const override;

private:
    /** Source particle container */
    WarpXParticleContainer* m_pc_src = nullptr;
    /** String containing species name of particles being written */
    std::string m_species_name;
};

#endif //
