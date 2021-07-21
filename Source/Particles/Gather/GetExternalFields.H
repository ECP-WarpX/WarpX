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
    const amrex::Real* AMREX_RESTRICT m_repeated_plasma_lens_starts = nullptr;
    const amrex::Real* AMREX_RESTRICT m_repeated_plasma_lens_lengths = nullptr;
    const amrex::Real* AMREX_RESTRICT m_repeated_plasma_lens_strengths = nullptr;
    int m_n_lenses;
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

            amrex::ParticleReal const uxp = m_ux[i];
            amrex::ParticleReal const uyp = m_uy[i];
            amrex::ParticleReal const uzp = m_uz[i];
            constexpr amrex::Real inv_c2 = 1._rt/(PhysConst::c*PhysConst::c);
            const amrex::Real inv_gamma = 1._rt/std::sqrt(1._rt + (uxp*uxp + uyp*uyp + uzp*uzp)*inv_c2);
            const amrex::ParticleReal vzp = uzp*inv_gamma;

            // This assumes that vzp > 0.
            amrex::ParticleReal const zl = z;
            amrex::ParticleReal const zr = z + vzp*m_dt;

            // This assumes that zl > 0.
            int i_lens = static_cast<int>(std::floor(zl/m_repeated_plasma_lens_period));
            i_lens = i_lens % m_n_lenses;
            amrex::Real const lens_start = m_repeated_plasma_lens_starts[i_lens] + i_lens*m_repeated_plasma_lens_period;
            amrex::Real const lens_end = lens_start + m_repeated_plasma_lens_lengths[i_lens];

            // Calculate the residence correction
            // frac will be 1 if the step is completely inside the lens, between 0 and 1
            // when entering or leaving the lens, and otherwise 0.
            amrex::Real fl = 0.;
            if (zl >= lens_start && zl < lens_end) fl = 1.;
            amrex::Real fr = 0.;
            if (zr >= lens_start && zr < lens_end) fr = 1.;
            amrex::Real frac = fl;
            amrex::Real dzi = 1./(vzp*m_dt);
            if (fl > fr) frac = (lens_end - zl)*dzi;
            if (fr > fl) frac = (zr - lens_start)*dzi;

            field_x += x*frac*m_repeated_plasma_lens_strengths[i_lens];
            field_y += y*frac*m_repeated_plasma_lens_strengths[i_lens];

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
