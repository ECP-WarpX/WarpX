#ifndef WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_
#define WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_

#include "Particles/WarpXParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <AMReX_REAL.H>

#include <limits>

enum ExternalFieldInitType { Constant, Parser };

/** \brief Base class for functors that assign external
 *         field values (E or B) to particles.
*/
struct GetExternalField
{
    ExternalFieldInitType m_type;

    amrex::GpuArray<amrex::ParticleReal, 3> m_field_value;

    ParserWrapper<4>* m_xfield_partparser = nullptr;
    ParserWrapper<4>* m_yfield_partparser = nullptr;
    ParserWrapper<4>* m_zfield_partparser = nullptr;
    GetParticlePosition m_get_position;
    amrex::Real m_time;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (long i,
                      amrex::ParticleReal& field_x,
                      amrex::ParticleReal& field_y,
                      amrex::ParticleReal& field_z) const noexcept
    {
        if (m_type == Constant)
        {
            field_x = m_field_value[0];
            field_y = m_field_value[1];
            field_z = m_field_value[2];
        }
        else if (m_type == Parser)
        {
            AMREX_ASSERT(m_xfield_partparser != nullptr);
            AMREX_ASSERT(m_yfield_partparser != nullptr);
            AMREX_ASSERT(m_zfield_partparser != nullptr);

            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);
            field_x = (*m_xfield_partparser)(x, y, z, m_time);
            field_y = (*m_yfield_partparser)(x, y, z, m_time);
            field_z = (*m_zfield_partparser)(x, y, z, m_time);
        }
        else
        {
            amrex::Abort("ExternalFieldInitType not known!!! \n");
        }
    }
};

/** \brief Functor that can be used to assign the external
 *         E field to a particle inside a ParallelFor kernel
*/
struct GetExternalEField : GetExternalField
{
    GetExternalEField () = default;

    GetExternalEField (const WarpXParIter& a_pti, int a_offset = 0) noexcept;
};

/** \brief Functor that can be used to assign the external
 *         B field to a particle inside a ParallelFor kernel
*/
struct GetExternalBField : GetExternalField
{
    GetExternalBField () = default;

    GetExternalBField (const WarpXParIter& a_pti, int a_offset = 0) noexcept;
};

#endif
