/* Copyright 2019 David Grote, Maxence Thevenet, Remi Lehe
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
#define WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_

#include "Particles/WarpXParticleContainer.H"

#include <AMReX.H>
#include <AMReX_REAL.H>

#include <limits>

/** \brief Extract the cartesian position coordinates of the particle
 *         p and store them in the variables `x`, `y`, `z`
 *         This is the implementaton and can't be called directly */
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void getParticlePosition (const P& p, const RType theta,
                          RType& x, RType& y, RType& z) const noexcept
{
#ifdef WARPX_DIM_RZ
    RType r = p.pos(0);
    x = r*std::cos(theta);
    y = r*std::sin(theta);
    z = p.pos(1);
#elif WARPX_DIM_3D
    x = p.pos(0);
    y = p.pos(1);
    z = p.pos(2);
#else
    x = p.pos(0);
    y = m_snan;
    z = p.pos(1);
#endif
}

/** \brief Functor that can be used to extract the positions of the macroparticles
 *         inside a ParallelFor kernel
 *
 * \param a_pti iterator to the tile containing the macroparticles
 * \param a_offset offset to apply to the particle indices
*/
struct GetParticlePosition
{
    using PType = WarpXParticleContainer::ParticleType;
    using SuperPType = WarpXParticleContainer::SuperParticleType;
    using RType = amrex::ParticleReal;

    const PType* AMREX_RESTRICT m_structs = nullptr;
#if (defined WARPX_DIM_RZ)
    const RType* m_theta = nullptr;
#elif (AMREX_SPACEDIM == 2)
    static constexpr RType m_snan = std::numeric_limits<RType>::quiet_NaN();
#endif

    GetParticlePosition () = default;

    GetParticlePosition (const WarpXParIter& a_pti, int a_offset = 0) noexcept
    {
        const auto& aos = a_pti.GetArrayOfStructs();
        m_structs = aos().dataPtr() + a_offset;
#if (defined WARPX_DIM_RZ)
        const auto& soa = a_pti.GetStructOfArrays();
        m_theta = soa.GetRealData(PIdx::theta).dataPtr() + a_offset;
#endif
    }

    /** \brief Extract the cartesian position coordinates of the particle
     *         located at index `i + a_offset` and store them in the variables
     *         `x`, `y`, `z` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const int i, RType& x, RType& y, RType& z) const noexcept
    {
        unpack_particle_position(m_structs[i], m_theta[i], x, y, z);
    }
};

/** \brief Functor that can be used to modify the positions of the macroparticles,
 *         inside a ParallelFor kernel.
 *
 * \param a_pti iterator to the tile being modified
 * \param a_offset offset to apply to the particle indices
*/
struct SetParticlePosition
{
    using PType = WarpXParticleContainer::ParticleType;
    using RType = amrex::ParticleReal;

    PType* AMREX_RESTRICT m_structs;
#if (defined WARPX_DIM_RZ)
    RType* AMREX_RESTRICT m_theta;
#endif
    SetParticlePosition (WarpXParIter& a_pti, int a_offset = 0) noexcept
    {
        auto& aos = a_pti.GetArrayOfStructs();
        m_structs = aos().dataPtr() + a_offset;
#if (defined WARPX_DIM_RZ)
        auto& soa = a_pti.GetStructOfArrays();
        m_theta = soa.GetRealData(PIdx::theta).dataPtr() + a_offset;
#endif
    }

    /** \brief Set the position of the particle at index `i + a_offset`
     *         to `x`, `y`, `z` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const int i, RType x, RType y, RType z) const noexcept
    {
#if defined(WARPX_DIM_XZ)
        amrex::ignore_unused(y);
#endif
#ifdef WARPX_DIM_RZ
        m_theta[i] = std::atan2(y, x);
        m_structs[i].pos(0) = std::sqrt(x*x + y*y);
        m_structs[i].pos(1) = z;
#elif WARPX_DIM_3D
        m_structs[i].pos(0) = x;
        m_structs[i].pos(1) = y;
        m_structs[i].pos(2) = z;
#else
        m_structs[i].pos(0) = x;
        m_structs[i].pos(1) = z;
#endif
    }
};

#endif // WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
