#ifndef WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_
#define WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_

#include "Particles/WarpXParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>

#include <limits>

using namespace amrex::literals;

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
    void operator () (long ip,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& field_Ez,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& field_Bz,
                      const amrex::Real gamma_boost,
                      const amrex::Real beta_boost,
                      const amrex::GpuArray<int, 3>& boost_direction) const noexcept
    {
        const amrex::Real c = PhysConst::c;
        const amrex::Real c2 = PhysConst::c * PhysConst::c;

        // E and B fields in the lab frame
        amrex::ParticleReal field_Exl = 0.0_rt, field_Eyl = 0.0_rt, field_Ezl = 0.0_rt;
        amrex::ParticleReal field_Bxl = 0.0_rt, field_Byl = 0.0_rt, field_Bzl = 0.0_rt;

        // E field in the lab frame
        if (m_E_type == Constant)
        {
            field_Exl = m_E_field_value[0];
            field_Eyl = m_E_field_value[1];
            field_Ezl = m_E_field_value[2];
        }
        else if (m_E_type == Parser)
        {
            AMREX_ASSERT(m_Ex_field_partparser != nullptr);
            AMREX_ASSERT(m_Ey_field_partparser != nullptr);
            AMREX_ASSERT(m_Ez_field_partparser != nullptr);

            // Space coordinates in the boosted frame
            amrex::ParticleReal xb, yb, zb;
            m_get_position(ip, xb, yb, zb);

            // Space-time coordinates in the lab frame
            amrex::Real tl;
            amrex::ParticleReal xl, yl, zl;

            // Projection of boosted frame space coordinates onto boost direction
            const amrex::Real boost_proj = xb * boost_direction[0] +
                                           yb * boost_direction[1] +
                                           zb * boost_direction[2];

            // Space-time coordinates in the lab frame
            tl = gamma_boost * (m_time + boost_proj * beta_boost / c);
            xl = xb + ((gamma_boost - 1.) * boost_proj +
                 gamma_boost * m_time * beta_boost * c) * boost_direction[0];
            yl = yb + ((gamma_boost - 1.) * boost_proj +
                 gamma_boost * m_time * beta_boost * c) * boost_direction[1];
            zl = zb + ((gamma_boost - 1.) * boost_proj +
                 gamma_boost * m_time * beta_boost * c) * boost_direction[2];

            field_Exl = (*m_Ex_field_partparser)(xl, yl, zl, tl);
            field_Eyl = (*m_Ey_field_partparser)(xl, yl, zl, tl);
            field_Ezl = (*m_Ez_field_partparser)(xl, yl, zl, tl);
        }
        else
        {
            amrex::Abort("ExternalFieldInitType not known!!! \n");
        }

        // B field in the lab frame
        if (m_B_type == Constant)
        {
            field_Bxl = m_B_field_value[0];
            field_Byl = m_B_field_value[1];
            field_Bzl = m_B_field_value[2];
        }
        else if (m_B_type == Parser)
        {
            AMREX_ASSERT(m_Bx_field_partparser != nullptr);
            AMREX_ASSERT(m_By_field_partparser != nullptr);
            AMREX_ASSERT(m_Bz_field_partparser != nullptr);

            // Space coordinates in the boosted frame
            amrex::ParticleReal xb, yb, zb;
            m_get_position(ip, xb, yb, zb);

            // Space-time coordinates in the lab frame
            amrex::Real tl;
            amrex::ParticleReal xl, yl, zl;

            // Projection of boosted frame space coordinates onto boost direction
            const amrex::Real boost_proj = xb * boost_direction[0] +
                                           yb * boost_direction[1] +
                                           zb * boost_direction[2];

            // Space-time coordinates in the lab frame
            tl = gamma_boost * (m_time + boost_proj * beta_boost / c);
            xl = xb + ((gamma_boost - 1.) * boost_proj +
                 gamma_boost * m_time * beta_boost * c) * boost_direction[0];
            yl = yb + ((gamma_boost - 1.) * boost_proj +
                 gamma_boost * m_time * beta_boost * c) * boost_direction[1];
            zl = zb + ((gamma_boost - 1.) * boost_proj +
                 gamma_boost * m_time * beta_boost * c) * boost_direction[2];

            field_Bxl = (*m_Bx_field_partparser)(xl, yl, zl, tl);
            field_Byl = (*m_By_field_partparser)(xl, yl, zl, tl);
            field_Bzl = (*m_Bz_field_partparser)(xl, yl, zl, tl);
        }
        else
        {
            amrex::Abort("ExternalFieldInitType not known!!! \n");
        }

        // Boost velocity
        const amrex::Vector<amrex::Real> v = {beta_boost * c * boost_direction[0],
                                              beta_boost * c * boost_direction[1],
                                              beta_boost * c * boost_direction[2]};

        const amrex::Vector<amrex::Real> v_cross_E = {v[1] * field_Ezl - v[2] * field_Eyl,
                                                      v[2] * field_Exl - v[0] * field_Ezl,
                                                      v[0] * field_Eyl - v[1] * field_Exl};

        const amrex::Vector<amrex::Real> v_cross_B = {v[1] * field_Bzl - v[2] * field_Byl,
                                                      v[2] * field_Bxl - v[0] * field_Bzl,
                                                      v[0] * field_Byl - v[1] * field_Bxl};

        const amrex::Real boost_dot_E = boost_direction[0] * field_Exl +
                                        boost_direction[1] * field_Eyl +
                                        boost_direction[2] * field_Ezl;

        const amrex::Real boost_dot_B = boost_direction[0] * field_Bxl +
                                        boost_direction[1] * field_Byl +
                                        boost_direction[2] * field_Bzl;

        // E field in the boosted frame
        field_Ex += gamma_boost * (field_Exl + v_cross_B[0]) - (gamma_boost - 1.) *
                    boost_dot_E * boost_direction[0];
        field_Ey += gamma_boost * (field_Eyl + v_cross_B[1]) - (gamma_boost - 1.) *
                    boost_dot_E * boost_direction[1];
        field_Ez += gamma_boost * (field_Ezl + v_cross_B[2]) - (gamma_boost - 1.) *
                    boost_dot_E * boost_direction[2];

        // B field in the boosted frame
        field_Bx += gamma_boost * (field_Bxl - v_cross_E[0] / c2) - (gamma_boost - 1.) *
                    boost_dot_B * boost_direction[0];
        field_By += gamma_boost * (field_Byl - v_cross_E[1] / c2) - (gamma_boost - 1.) *
                    boost_dot_B * boost_direction[1];
        field_Bz += gamma_boost * (field_Bzl - v_cross_E[2] / c2) - (gamma_boost - 1.) *
                    boost_dot_B * boost_direction[2];
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
