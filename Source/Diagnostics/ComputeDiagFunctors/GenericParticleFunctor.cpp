/* Copyright 2023 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "GenericParticleFunctor.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BaseFwd.H>

/**
 * \brief Functor to select particles to be written and copy them to `pc_src`
 */
GenericParticleFunctor::GenericParticleFunctor (
                              WarpXParticleContainer *pc_src,
                              std::string species_name )
    : m_pc_src{pc_src}, m_species_name{std::move(species_name)}
{

}

void
GenericParticleFunctor::operator () (PinnedMemoryParticleContainer& /*pc_dst*/, int& /*totalParticleCounter*/, int /*i_buffer*/) const
{
}
