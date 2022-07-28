#ifndef WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_
#define WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_

#include "Particles/Pusher/GetAndSetPosition.H"

#include "Particles/WarpXParticleContainer_fwd.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

enum ExternalFieldInitType { None, Constant, Parser, RepeatedPlasmaLens, Unknown };

/** \brief Functor class that assigns external
 *         field values (E and B) to particles.
*/
struct GetExternalEBField
{

    GetExternalEBField () = default;

    GetExternalEBField (const WarpXParIter& a_pti, int a_offset = 0) noexcept;

    ExternalFieldInitType m_Etype;
    ExternalFieldInitType m_Btype;

    amrex::ParticleReal m_gamma_boost;
    amrex::ParticleReal m_uz_boost;

    amrex::GpuArray<amrex::ParticleReal, 3> m_Efield_value;
    amrex::GpuArray<amrex::ParticleReal, 3> m_Bfield_value;

    amrex::ParserExecutor<4> m_Exfield_partparser;
    amrex::ParserExecutor<4> m_Eyfield_partparser;
    amrex::ParserExecutor<4> m_Ezfield_partparser;
    amrex::ParserExecutor<4> m_Bxfield_partparser;
    amrex::ParserExecutor<4> m_Byfield_partparser;
    amrex::ParserExecutor<4> m_Bzfield_partparser;

    GetParticlePosition m_get_position;
    amrex::Real m_time;

    amrex::ParticleReal m_repeated_plasma_lens_period;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_starts = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_lengths = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_strengths_E = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_strengths_B = nullptr;
    int m_n_lenses;
    amrex::Real m_dt;
    const amrex::ParticleReal* AMREX_RESTRICT m_ux = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uy = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uz = nullptr;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (long i,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& field_Ez,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& field_Bz) const noexcept
    {
        using namespace amrex::literals;

        if (m_Etype == None && m_Btype == None) return;

        amrex::ParticleReal Ex = 0._prt;
        amrex::ParticleReal Ey = 0._prt;
        amrex::ParticleReal Ez = 0._prt;
        amrex::ParticleReal Bx = 0._prt;
        amrex::ParticleReal By = 0._prt;
        amrex::ParticleReal Bz = 0._prt;

        constexpr amrex::ParticleReal inv_c2 = 1._prt/(PhysConst::c*PhysConst::c);

        if (m_Etype == Constant)
        {
            Ex = m_Efield_value[0];
            Ey = m_Efield_value[1];
            Ez = m_Efield_value[2];
        }
        else if (m_Etype == ExternalFieldInitType::Parser)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);
            amrex::Real lab_time = m_time;
            if (m_gamma_boost > 1._prt) {
                lab_time = m_gamma_boost*m_time + m_uz_boost*z*inv_c2;
                z = m_gamma_boost*z + m_uz_boost*m_time;
            }
            Ex = m_Exfield_partparser((amrex::ParticleReal) x, (amrex::ParticleReal) y, (amrex::ParticleReal) z, lab_time);
            Ey = m_Eyfield_partparser((amrex::ParticleReal) x, (amrex::ParticleReal) y, (amrex::ParticleReal) z, lab_time);
            Ez = m_Ezfield_partparser((amrex::ParticleReal) x, (amrex::ParticleReal) y, (amrex::ParticleReal) z, lab_time);
        }

        if (m_Btype == Constant)
        {
            Bx = m_Bfield_value[0];
            By = m_Bfield_value[1];
            Bz = m_Bfield_value[2];
        }
        else if (m_Btype == ExternalFieldInitType::Parser)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);
            amrex::Real lab_time = m_time;
            if (m_gamma_boost > 1._prt) {
                lab_time = m_gamma_boost*m_time + m_uz_boost*z*inv_c2;
                z = m_gamma_boost*z + m_uz_boost*m_time;
            }
            Bx = m_Bxfield_partparser(x, y, z, lab_time);
            By = m_Byfield_partparser(x, y, z, lab_time);
            Bz = m_Bzfield_partparser(x, y, z, lab_time);
        }

        if (m_Etype == RepeatedPlasmaLens ||
            m_Btype == RepeatedPlasmaLens)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);

            const amrex::ParticleReal uxp = m_ux[i];
            const amrex::ParticleReal uyp = m_uy[i];
            const amrex::ParticleReal uzp = m_uz[i];

            const amrex::ParticleReal gamma = std::sqrt(1._prt + (uxp*uxp + uyp*uyp + uzp*uzp)*inv_c2);
            const amrex::ParticleReal vzp = uzp/gamma;

            amrex::ParticleReal zl = z;
            amrex::ParticleReal zr = z + vzp*m_dt;

            if (m_gamma_boost > 1._prt) {
                zl = m_gamma_boost*zl + m_uz_boost*m_time;
                zr = m_gamma_boost*zr + m_uz_boost*(m_time + m_dt);
            }

            // This assumes that zl > 0.
            int i_lens = static_cast<int>(std::floor(zl/m_repeated_plasma_lens_period));
            i_lens = i_lens % m_n_lenses;
            amrex::ParticleReal const lens_start = m_repeated_plasma_lens_starts[i_lens] + i_lens*m_repeated_plasma_lens_period;
            amrex::ParticleReal const lens_end = lens_start + m_repeated_plasma_lens_lengths[i_lens];

            // Calculate the residence correction
            // frac will be 1 if the step is completely inside the lens, between 0 and 1
            // when entering or leaving the lens, and otherwise 0.
            // This assumes that vzp > 0.
            amrex::ParticleReal fl = 0.;
            if (zl >= lens_start && zl < lens_end) fl = 1.;
            amrex::ParticleReal fr = 0.;
            if (zr >= lens_start && zr < lens_end) fr = 1.;
            amrex::ParticleReal frac = fl;
            if (fl > fr) frac = (lens_end - zl)/(zr - zl);
            if (fr > fl) frac = (zr - lens_start)/(zr - zl);

            // Note that "+=" is used since the fields may have been set above
            // if a different E or Btype was specified.
            Ex += x*frac*m_repeated_plasma_lens_strengths_E[i_lens];
            Ey += y*frac*m_repeated_plasma_lens_strengths_E[i_lens];
            Bx += +y*frac*m_repeated_plasma_lens_strengths_B[i_lens];
            By += -x*frac*m_repeated_plasma_lens_strengths_B[i_lens];

        }

        if (m_gamma_boost > 1._prt) {
            // Transform the fields to the boosted frame
            const amrex::ParticleReal Ex_boost = m_gamma_boost*Ex - m_uz_boost*By;
            const amrex::ParticleReal Ey_boost = m_gamma_boost*Ey + m_uz_boost*Bx;
            const amrex::ParticleReal Bx_boost = m_gamma_boost*Bx + m_uz_boost*Ey*inv_c2;
            const amrex::ParticleReal By_boost = m_gamma_boost*By - m_uz_boost*Ex*inv_c2;
            Ex = Ex_boost;
            Ey = Ey_boost;
            Bx = Bx_boost;
            By = By_boost;
        }

        field_Ex += Ex;
        field_Ey += Ey;
        field_Ez += Ez;
        field_Bx += Bx;
        field_By += By;
        field_Bz += Bz;

    }
};

#endif
