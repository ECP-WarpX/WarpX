#ifndef FILTERFUNCTORS_H
#define FILTERFUNCTORS_H

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "WarpXParticleContainer.H"

using ParticleType = typename WarpXParticleContainer::ParticleType;

/**
 * \brief Functor that returns 0 or 1 depending on a random draw per particle
 */
struct RandomFilter
{
public:
    /** constructor
     * \param a_is_active whether the test is active
     * \param a_fraction fraction of particles to select
     */
    RandomFilter(bool a_is_active, amrex::Real a_fraction)
        : m_is_active(a_is_active), m_fraction(a_fraction) {};

    /**
     * \brief draw random number, return 1 if number < m_fraction, 1 otherwise
     * \param p one particle
     * \return whether or not the particle is selected
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator () (const ParticleType& p) const noexcept
    {
        if ( !m_is_active )
            return 1;
        if ( amrex::Random() < m_fraction )
            return 1;
        return 0;
    }
private:
    bool m_is_active;
    amrex::Real m_fraction;
};

/**
 * \brief Functor that returns 1 if stride divide particle_id, 0 otherwise
 */
struct UniformFilter
{
public:
    /** constructor
     * \param a_is_active whether the test is active
     * \param a_stride one particle every a_stride is written to file
     */
    UniformFilter(bool a_is_active, long a_stride)
        : m_is_active(a_is_active), m_stride(a_stride) {};

    /**
     * \brief return 1 if stride divide particle_id, 0 otherwise
     * \param p one particle
     * \return whether or not the particle is selected
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator () (const ParticleType& p) const noexcept
    {
        if ( !m_is_active )
            return 1;
        if ( p.id()%m_stride == 0 )
            return 1;
        return 0;
    }
private:
    bool m_is_active;
    long m_stride;
};

#endif // FILTERFUNCTORS_H
