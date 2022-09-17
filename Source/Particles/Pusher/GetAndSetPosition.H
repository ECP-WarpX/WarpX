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

#include <cmath>
#include <limits>

/** \brief Extract the cartesian position coordinates of the particle
 *         p and store them in the variables `x`, `y`, `z`
 *         This version operates on a SuperParticle */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_particle_position (const WarpXParticleContainer::SuperParticleType& p,
                            amrex::ParticleReal& x,
                            amrex::ParticleReal& y,
                            amrex::ParticleReal& z) noexcept
{
#ifdef WARPX_DIM_RZ
    amrex::ParticleReal theta = p.rdata(PIdx::theta);
    amrex::ParticleReal r = p.pos(0);
    x = r*std::cos(theta);
    y = r*std::sin(theta);
    z = p.pos(1);
#elif WARPX_DIM_3D
    x = p.pos(0);
    y = p.pos(1);
    z = p.pos(2);
#elif WARPX_DIM_XZ
    x = p.pos(0);
    y = amrex::ParticleReal(0.0);
    z = p.pos(1);
#else
    x = amrex::ParticleReal(0.0);
    y = amrex::ParticleReal(0.0);
    z = p.pos(0);
#endif
}

/** \brief Functor that can be used to extract the positions of the macroparticles
 *         inside a ParallelFor kernel
*/
struct GetParticlePosition
{
    using RType = amrex::ParticleReal;

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
    RType* AMREX_RESTRICT m_x = nullptr;
    RType* AMREX_RESTRICT m_z = nullptr;
#elif defined(WARPX_DIM_3D)
    const RType* AMREX_RESTRICT m_x = nullptr;
    const RType* AMREX_RESTRICT m_y = nullptr;
    const RType* AMREX_RESTRICT m_z = nullptr;
#elif defined(WARPX_DIM_1D_Z)
    RType* AMREX_RESTRICT m_z = nullptr;;
#endif

#if defined(WARPX_DIM_RZ)
    const RType* m_theta = nullptr;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    static constexpr RType m_y_default = RType(0.0);
#elif defined(WARPX_DIM_1D_Z)
    static constexpr RType m_x_default = RType(0.0);
    static constexpr RType m_y_default = RType(0.0);
#endif

    GetParticlePosition () = default;

    /** Constructor
     *
     * \tparam ptiType the type of the particle iterator used in the constructor
     *
     * \param a_pti iterator to the tile containing the macroparticles
     * \param a_offset offset to apply to the particle indices
     */
    template <typename ptiType>
    GetParticlePosition (const ptiType& a_pti, int a_offset = 0) noexcept
    {
        const auto& soa = a_pti.GetStructOfArrays();

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
        m_x = soa.GetRealData(0).dataPtr() + a_offset;
        m_z = soa.GetRealData(1).dataPtr() + a_offset;
#elif defined(WARPX_DIM_3D)
        m_x = soa.GetRealData(0).dataPtr() + a_offset;
        m_y = soa.GetRealData(1).dataPtr() + a_offset;
        m_z = soa.GetRealData(2).dataPtr() + a_offset;
#elif defined(WARPX_DIM_1D_Z)
        m_z = soa.GetRealData(0).dataPtr() + a_offset;
#endif
#if defined(WARPX_DIM_RZ)
        m_theta = soa.GetRealData(PIdx::theta).dataPtr() + a_offset;
#endif
    }

    /** \brief Extract the cartesian position coordinates of the particle
     *         located at index `i + a_offset` and store them in the variables
     *         `x`, `y`, `z` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const int i, RType& x, RType& y, RType& z) const noexcept
    {
#ifdef WARPX_DIM_RZ
        RType r = m_x[i];
        x = r*std::cos(m_theta[i]);
        y = r*std::sin(m_theta[i]);
        z = m_z[i];
#elif WARPX_DIM_3D
        x = m_x[i];
        y = m_y[i];
        z = m_z[i];
#elif WARPX_DIM_XZ
        x = m_x[i];
        y = m_y_default;
        z = m_z[i];
#else
        x = m_x_default;
        y = m_y_default;
        z = m_z[i];
#endif
    }

    /** \brief Extract the position coordinates of the particle as stored
     *         located at index `i + a_offset` and store them in the variables
     *         `x`, `y`, `z`
     *         This is only different for RZ since this returns (r, theta, z)
     *         in that case. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void AsStored (const int i, RType& x, RType& y, RType& z) const noexcept
    {
#ifdef WARPX_DIM_RZ
        x = m_x[i];
        y = m_theta[i];
        z = m_z[i];
#elif WARPX_DIM_3D
        x = m_x[i];
        y = m_y[i];
        z = m_z[i];
#elif WARPX_DIM_XZ
        x = m_x[i];
        y = m_y_default;
        z = m_z[i];
#else
        x = m_x_default;
        y = m_y_default;
        z = m_z[i];
#endif
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
    using RType = amrex::ParticleReal;

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
    RType* AMREX_RESTRICT m_x;
    RType* AMREX_RESTRICT m_z;
#elif defined(WARPX_DIM_3D)
    RType* AMREX_RESTRICT m_x;
    RType* AMREX_RESTRICT m_y;
    RType* AMREX_RESTRICT m_z;
#elif defined(WARPX_DIM_1D_Z)
    RType* AMREX_RESTRICT m_z;
#endif
#if defined(WARPX_DIM_RZ)
    RType* AMREX_RESTRICT m_theta;
#endif

    template <typename ptiType>
    SetParticlePosition (const ptiType& a_pti, int a_offset = 0) noexcept
    {
        auto& soa = a_pti.GetStructOfArrays();
#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
        m_x = soa.GetRealData(0).dataPtr() + a_offset;
        m_z = soa.GetRealData(1).dataPtr() + a_offset;
#elif defined(WARPX_DIM_3D)
        m_x = soa.GetRealData(0).dataPtr() + a_offset;
        m_y = soa.GetRealData(1).dataPtr() + a_offset;
        m_z = soa.GetRealData(2).dataPtr() + a_offset;
#elif defined(WARPX_DIM_1D_Z)
        m_z = soa.GetRealData(0).dataPtr() + a_offset;
#endif
#if defined(WARPX_DIM_RZ)
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
#if defined(WARPX_DIM_1D_Z)
        amrex::ignore_unused(x,y);
#endif
#ifdef WARPX_DIM_RZ
        m_theta[i] = std::atan2(y, x);
        m_x[i] = std::sqrt(x*x + y*y);
        m_z[i] = z;
#elif WARPX_DIM_3D
        m_x[i] = x;
        m_y[i] = y;
        m_z[i] = z;
#elif WARPX_DIM_XZ
        m_x[i] = x;
        m_z[i] = z;
#else
        m_z[i] = z;
#endif
    }

    /** \brief Set the position of the particle at index `i + a_offset`
     *         to `x`, `y`, `z`
     *         This is only different for RZ since the input should
     *         be (r, theta, z) in that case. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void AsStored (const int i, RType x, RType y, RType z) const noexcept
    {
#if defined(WARPX_DIM_XZ)
        amrex::ignore_unused(y);
#endif
#if defined(WARPX_DIM_1D_Z)
        amrex::ignore_unused(x,y);
#endif
#ifdef WARPX_DIM_RZ
        m_x[i] = x;
        m_theta[i] = y;
        m_z[i] = z;
#elif WARPX_DIM_3D
        m_x[i] = x;
        m_y[i] = y;
        m_z[i] = z;
#elif WARPX_DIM_XZ
        m_x[i] = x;
        m_z[i] = z;
#else
        m_z[i] = z;
#endif
    }
};

#endif // WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
