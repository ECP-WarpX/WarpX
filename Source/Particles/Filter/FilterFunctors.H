#ifndef FILTERFUNCTORS_H
#define FILTERFUNCTORS_H

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "WarpXParticleContainer.H"
#include "WarpX.H"

using SuperParticleType = typename WarpXParticleContainer::SuperParticleType;

/**
 * \brief Functor that returns 0 or 1 depending on a random draw per particle
 */
struct RandomFilter
{
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
    bool operator () (const SuperParticleType& p) const noexcept
    {
        if ( !m_is_active )
            return 1;
        if ( amrex::Random() < m_fraction )
            return 1;
        return 0;
    }
private:
    bool m_is_active = false; //! select all particles if false
    amrex::Real m_fraction = 1.0; //! range: [0.0:1.0] where 0 is no & 1 is all particles
};

/**
 * \brief Functor that returns 1 if stride divide particle_id, 0 otherwise
 */
struct UniformFilter
{
    /** constructor
     * \param a_is_active whether the test is active
     * \param a_stride one particle every a_stride is written to file
     */
    UniformFilter(bool a_is_active, int a_stride)
        : m_is_active(a_is_active), m_stride(a_stride) {};

    /**
     * \brief return 1 if stride divide particle_id, 0 otherwise
     * \param p one particle
     * \return whether or not the particle is selected
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator () (const SuperParticleType& p) const noexcept
    {
        if ( !m_is_active )
            return 1;
        if ( p.id()%m_stride == 0 )
            return 1;
        return 0;
    }
private:
    bool m_is_active = false; //! select all particles if false
    int m_stride = 0; //! selection of every n-th particle
};

/**
 * \brief Functor that returns 0 or 1 depending on a parser selection
 */
struct ParserFilter
{
    /** constructor
     * \param a_is_active whether the test is active
     */
    ParserFilter(bool a_is_active, int a_species_index)
        : m_is_active(a_is_active), m_species_index(a_species_index){};

    /**
     * \brief return 1 if the particle is selected by the parser
     * \param p one particle
     * \return whether or not the particle is selected
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator () (const SuperParticleType& p) const noexcept
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(m_species_index);
        ParserWrapper<7> *function_partparser =
            myspc.m_particle_filter_parser.get();
        amrex::Real t = WarpX::GetInstance().gett_new(0);
        amrex::Real x  = p.pos(0);
        amrex::Real y  = p.pos(1);
        amrex::Real z  = p.pos(2);
        amrex::Real ux = p.rdata(PIdx::ux);
        amrex::Real uy = p.rdata(PIdx::uy);
        amrex::Real uz = p.rdata(PIdx::uz);
        if ( !m_is_active )
            return 1;
        if ( (*function_partparser)(t,x,y,z,ux,uy,uz) != 0.0 )
            return 1;
        return 0;
    }
private:
    bool m_is_active = false; //! select all particles if false
    int m_species_index = 0;  //! species index
};

#endif // FILTERFUNCTORS_H
