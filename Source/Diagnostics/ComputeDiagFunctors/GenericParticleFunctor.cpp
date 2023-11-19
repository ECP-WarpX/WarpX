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

SelectParticles::SelectParticles (const WarpXParIter& /*a_pti*/, TmpParticles& /*tmp_particle_data*/,
                                  amrex::Real /*current_z_boost*/, amrex::Real /*old_z_boost*/,
                                  int /*a_offset*/)
{
//    m_get_position = GetParticlePosition<PIdx>(a_pti, a_offset);

//    const auto lev = a_pti.GetLevel();
//    const auto index = a_pti.GetPairIndex();
}

/**
 * \brief Functor to compute Lorentz Transform and store the selected particles in existing
 * particle buffers
 */
GenericParticleFunctor::GenericParticleFunctor (
                              WarpXParticleContainer *pc_src,
                              std::string species_name,
                              int /*num_buffers*/)
    : m_pc_src{pc_src}, m_species_name{std::move(species_name)}
{
    InitData();
}

void
GenericParticleFunctor::operator () (PinnedMemoryParticleContainer& /*pc_dst*/, int& /*totalParticleCounter*/, int /*i_buffer*/) const
{
}

void
GenericParticleFunctor::PrepareFunctorData ( int /*i_buffer*/, bool /*z_slice_in_domain*/,
                              amrex::Real /*old_z_boost*/, amrex::Real /*current_z_boost*/,
                              amrex::Real /*t_lab*/, int /*snapshot_full*/)
{
}
