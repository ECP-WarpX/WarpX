/* Copyright 2020 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_COPYPARTICLEATTRIBS_H_
#define WARPX_PARTICLES_PUSHER_COPYPARTICLEATTRIBS_H_

#include "Particles/WarpXParticleContainer.H"

#include <AMReX_REAL.H>

#include <limits>


/** \brief Functor that creates copies of the current particle
 *         positions and momenta for later use. This is needed
 *         by the back-transformed diagnostics.
 *
 * \param a_pti iterator to the tile containing the macroparticles
 * \param a_tmp holder for the temporary particle data
*/
struct CopyParticleAttribs
{
    using TmpParticles = WarpXParticleContainer::TmpParticles;

    GetParticlePosition m_get_position;

    const amrex::ParticleReal* AMREX_RESTRICT uxp = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT uyp = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT uzp = nullptr;

    amrex::ParticleReal* AMREX_RESTRICT xpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT ypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT zpold = nullptr;

    amrex::ParticleReal* AMREX_RESTRICT uxpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT uypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT uzpold = nullptr;

    CopyParticleAttribs (const WarpXParIter& a_pti, TmpParticles& tmp_particle_data) noexcept
    {
        if (tmp_particle_data.size() == 0) return;

        auto& attribs = a_pti.GetAttribs();
        uxp = attribs[PIdx::ux].dataPtr();
        uyp = attribs[PIdx::uy].dataPtr();
        uzp = attribs[PIdx::uz].dataPtr();

        const auto lev = a_pti.GetLevel();
        const auto index = a_pti.GetPairIndex();
        xpold  = tmp_particle_data[lev][index][TmpIdx::xold ].dataPtr();
        ypold  = tmp_particle_data[lev][index][TmpIdx::yold ].dataPtr();
        zpold  = tmp_particle_data[lev][index][TmpIdx::zold ].dataPtr();
        uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
        uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
        uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();

        m_get_position = GetParticlePosition(a_pti);
    }

    /** \brief copy the position and momentum of particle i to the
     *         temporary data holder
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const int i) const noexcept
    {
        AMREX_ASSERT(uxp != nullptr);
        AMREX_ASSERT(uyp != nullptr);
        AMREX_ASSERT(uzp != nullptr);

        AMREX_ASSERT(xpold != nullptr);
        AMREX_ASSERT(ypold != nullptr);
        AMREX_ASSERT(zpold != nullptr);

        AMREX_ASSERT(uxpold != nullptr);
        AMREX_ASSERT(uypold != nullptr);
        AMREX_ASSERT(uzpold != nullptr);

        amrex::ParticleReal x, y, z;
        m_get_position(i, x, y, z);

        xpold[i] = x;
        ypold[i] = y;
        zpold[i] = z;

        uxpold[i] = uxp[i];
        uypold[i] = uyp[i];
        uzpold[i] = uzp[i];
    }
};

#endif // WARPX_PARTICLES_PUSHER_COPYPARTICLEATTRIBS_H_
