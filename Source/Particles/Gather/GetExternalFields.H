#ifndef WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_
#define WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_

#include "Particles/WarpXParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <AMReX_REAL.H>

#include <limits>

enum ExternalFieldInitType { Constant, Parser };

/** \brief Base class for functors that assign the external E and B fields
 *         to a particle inside a ParallelFor kernel
 */
struct GetExternalField
{
    ExternalFieldInitType m_E_type;
    ExternalFieldInitType m_B_type;

    amrex::GpuArray<amrex::ParticleReal, 3> m_E_field_value;
    amrex::GpuArray<amrex::ParticleReal, 3> m_B_field_value;

    ParserWrapper<4>* m_Ex_field_partparser = nullptr;
    ParserWrapper<4>* m_Ey_field_partparser = nullptr;
    ParserWrapper<4>* m_Ez_field_partparser = nullptr;

    ParserWrapper<4>* m_Bx_field_partparser = nullptr;
    ParserWrapper<4>* m_By_field_partparser = nullptr;
    ParserWrapper<4>* m_Bz_field_partparser = nullptr;

    GetParticlePosition m_get_position;
    amrex::Real m_time;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (long i,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& field_Ez,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& field_Bz,
                      amrex::Real gamma_boost,
                      amrex::Real beta_boost) const noexcept
    {
        // Lab frame
        if (gamma_boost == 1.)
        {
            // E field
            if (m_E_type == Constant)
            {
                field_Ex += m_E_field_value[0];
                field_Ey += m_E_field_value[1];
                field_Ez += m_E_field_value[2];
            }
            else if (m_E_type == Parser)
            {
                AMREX_ASSERT(m_Ex_field_partparser != nullptr);
                AMREX_ASSERT(m_Ey_field_partparser != nullptr);
                AMREX_ASSERT(m_Ez_field_partparser != nullptr);

                amrex::ParticleReal x, y, z;
                m_get_position(i, x, y, z);
                field_Ex += (*m_Ex_field_partparser)(x, y, z, m_time);
                field_Ey += (*m_Ey_field_partparser)(x, y, z, m_time);
                field_Ez += (*m_Ez_field_partparser)(x, y, z, m_time);
            }
            else
            {
                amrex::Abort("ExternalFieldInitType not known!!! \n");
            }

            // B field
            if (m_B_type == Constant)
            {
                field_Bx += m_B_field_value[0];
                field_By += m_B_field_value[1];
                field_Bz += m_B_field_value[2];
            }
            else if (m_B_type == Parser)
            {
                AMREX_ASSERT(m_Bx_field_partparser != nullptr);
                AMREX_ASSERT(m_By_field_partparser != nullptr);
                AMREX_ASSERT(m_Bz_field_partparser != nullptr);

                amrex::ParticleReal x, y, z;
                m_get_position(i, x, y, z);
                field_Bx += (*m_Bx_field_partparser)(x, y, z, m_time);
                field_By += (*m_By_field_partparser)(x, y, z, m_time);
                field_Bz += (*m_Bz_field_partparser)(x, y, z, m_time);
            }
            else
            {
                amrex::Abort("ExternalFieldInitType not known!!! \n");
            }
        }

        // Boosted frame
        else if (gamma_boost > 1.)
        {
        }
    }
};

/** \brief Functor that can be used to assign the external E and B fields
 *         to a particle inside a ParallelFor kernel
 */
struct GetExternalEBField : GetExternalField
{
    GetExternalEBField () = default;

    GetExternalEBField (const WarpXParIter& a_pti, int a_offset = 0) noexcept;
};

#endif
