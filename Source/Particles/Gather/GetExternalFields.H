#ifndef WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_
#define WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_

#include "Particles/WarpXParticleContainer.H"

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
    GetExternalEField (const WarpXParIter& a_pti, int a_offset = 0) noexcept
    {
        auto& warpx = WarpX::GetInstance();
        auto& mypc = warpx.GetPartContainer();
        if (mypc.m_E_ext_particle_s=="constant" || mypc.m_E_ext_particle_s=="default")
        {
            m_type = Constant;
            m_field_value[0] = mypc.m_E_external_particle[0];
            m_field_value[1] = mypc.m_E_external_particle[1];
            m_field_value[2] = mypc.m_E_external_particle[2];
        }
        else if (mypc.m_E_ext_particle_s=="parse_e_ext_particle_function")
        {
            m_type = Parser;
            m_time = warpx.gett_new(a_pti.GetLevel());
            m_get_position = GetParticlePosition(a_pti, a_offset);
            m_xfield_partparser = mypc.m_Ex_particle_parser.get();
            m_yfield_partparser = mypc.m_Ey_particle_parser.get();
            m_zfield_partparser = mypc.m_Ez_particle_parser.get();
        }
    }
};

/** \brief Functor that can be used to assign the external
 *         B field to a particle inside a ParallelFor kernel
*/
struct GetExternalBField : GetExternalField
{
    GetExternalBField (const WarpXParIter& a_pti, int a_offset = 0) noexcept
    {
        auto& warpx = WarpX::GetInstance();
        auto& mypc = warpx.GetPartContainer();
        if (mypc.m_B_ext_particle_s=="constant" || mypc.m_B_ext_particle_s=="default")
        {
            m_type = Constant;
            m_field_value[0] = mypc.m_B_external_particle[0];
            m_field_value[1] = mypc.m_B_external_particle[1];
            m_field_value[2] = mypc.m_B_external_particle[2];
        }
        else if (mypc.m_B_ext_particle_s=="parse_e_ext_particle_function")
        {
            m_type = Parser;
            m_time = warpx.gett_new(a_pti.GetLevel());
            m_get_position = GetParticlePosition(a_pti, a_offset);
            m_xfield_partparser = mypc.m_Bx_particle_parser.get();
            m_yfield_partparser = mypc.m_By_particle_parser.get();
            m_zfield_partparser = mypc.m_Bz_particle_parser.get();
        }
    }
};

#endif
