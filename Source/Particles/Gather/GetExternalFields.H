#ifndef WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_
#define WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_

#include "Particles/Pusher/GetAndSetPosition.H"

#include "Particles/WarpXParticleContainer_fwd.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

enum ExternalFieldInitType { Constant, Parser, RepeatedPlasmaLens };

/** \brief Base class for functors that assign external
 *         field values (E or B) to particles.
*/
struct GetExternalField
{
    ExternalFieldInitType m_type;

    amrex::GpuArray<amrex::ParticleReal, 3> m_field_value;

    amrex::ParserExecutor<4> m_xfield_partparser;
    amrex::ParserExecutor<4> m_yfield_partparser;
    amrex::ParserExecutor<4> m_zfield_partparser;
    GetParticlePosition m_get_position;
    amrex::Real m_time;

    amrex::Real m_repeated_plasma_lens_period;
    amrex::Real m_gamma_boost;
    amrex::Real m_uz_boost;
    const amrex::Real* AMREX_RESTRICT m_repeated_plasma_lens_starts = nullptr;
    const amrex::Real* AMREX_RESTRICT m_repeated_plasma_lens_lengths = nullptr;
    const amrex::Real* AMREX_RESTRICT m_repeated_plasma_lens_strengths_E = nullptr;
    const amrex::Real* AMREX_RESTRICT m_repeated_plasma_lens_strengths_B = nullptr;
    int m_n_lenses;
    int m_lens_is_electric;
    amrex::Real m_dt;
    const amrex::ParticleReal* AMREX_RESTRICT m_ux = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uy = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uz = nullptr;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (long i,
                      amrex::ParticleReal& field_x,
                      amrex::ParticleReal& field_y,
                      amrex::ParticleReal& field_z) const noexcept
    {
        using namespace amrex::literals;
        if (m_type == Constant)
        {
            field_x += m_field_value[0];
            field_y += m_field_value[1];
            field_z += m_field_value[2];
        }
        else if (m_type == Parser)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);
            field_x += m_xfield_partparser(x, y, z, m_time);
            field_y += m_yfield_partparser(x, y, z, m_time);
            field_z += m_zfield_partparser(x, y, z, m_time);
        }
        else if (m_type == RepeatedPlasmaLens)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);

            const amrex::ParticleReal uxp = m_ux[i];
            const amrex::ParticleReal uyp = m_uy[i];
            const amrex::ParticleReal uzp = m_uz[i];

            constexpr amrex::Real inv_c2 = 1._rt/(PhysConst::c*PhysConst::c);
            const amrex::ParticleReal gamma = std::sqrt(1._rt + (uxp*uxp + uyp*uyp + uzp*uzp)*inv_c2);
            const amrex::ParticleReal vzp = uzp/gamma;

            amrex::ParticleReal zl = z;
            amrex::ParticleReal zr = z + vzp*m_dt;

            if (m_gamma_boost > 1._rt) {
                zl = m_gamma_boost*zl + m_uz_boost*m_time;
                zr = m_gamma_boost*zr + m_uz_boost*m_time;
            }

            // This assumes that zl > 0.
            int i_lens = static_cast<int>(std::floor(zl/m_repeated_plasma_lens_period));
            i_lens = i_lens % m_n_lenses;
            amrex::Real const lens_start = m_repeated_plasma_lens_starts[i_lens] + i_lens*m_repeated_plasma_lens_period;
            amrex::Real const lens_end = lens_start + m_repeated_plasma_lens_lengths[i_lens];

            // Calculate the residence correction
            // frac will be 1 if the step is completely inside the lens, between 0 and 1
            // when entering or leaving the lens, and otherwise 0.
            // This assumes that vzp > 0.
            amrex::Real fl = 0.;
            if (zl >= lens_start && zl < lens_end) fl = 1.;
            amrex::Real fr = 0.;
            if (zr >= lens_start && zr < lens_end) fr = 1.;
            amrex::Real frac = fl;
            if (fl > fr) frac = (lens_end - zl)/(zr - zl);
            if (fr > fl) frac = (zr - lens_start)/(zr - zl);

            if (m_lens_is_electric) {
                amrex::Real Ex = x*frac*m_repeated_plasma_lens_strengths_E[i_lens];
                amrex::Real Ey = y*frac*m_repeated_plasma_lens_strengths_E[i_lens];
                if (m_gamma_boost > 1._rt) {
                    // Transform the fields to the boosted frame
                    const amrex::Real Bx = +y*frac*m_repeated_plasma_lens_strengths_B[i_lens];
                    const amrex::Real By = -x*frac*m_repeated_plasma_lens_strengths_B[i_lens];
                    const amrex::Real vz_boost = m_uz_boost/m_gamma_boost;
                    const amrex::Real Ex_boost = m_gamma_boost*(Ex - vz_boost*By);
                    const amrex::Real Ey_boost = m_gamma_boost*(Ey + vz_boost*Bx);
                    Ex = Ex_boost;
                    Ey = Ey_boost;
                }
                field_x += Ex;
                field_y += Ey;
            } else {
                amrex::Real Bx = +y*frac*m_repeated_plasma_lens_strengths_B[i_lens];
                amrex::Real By = -x*frac*m_repeated_plasma_lens_strengths_B[i_lens];
                if (m_gamma_boost > 1._rt) {
                    // Transform the fields to the boosted frame
                    const amrex::Real Ex = x*frac*m_repeated_plasma_lens_strengths_E[i_lens];
                    const amrex::Real Ey = y*frac*m_repeated_plasma_lens_strengths_E[i_lens];
                    const amrex::Real vz_boost = m_uz_boost/m_gamma_boost;
                    const amrex::Real Bx_boost = m_gamma_boost*(Bx + vz_boost*Ey*inv_c2);
                    const amrex::Real By_boost = m_gamma_boost*(By - vz_boost*Ex*inv_c2);
                    Bx = Bx_boost;
                    By = By_boost;
                }
                field_x += Bx;
                field_y += By;
            }

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
