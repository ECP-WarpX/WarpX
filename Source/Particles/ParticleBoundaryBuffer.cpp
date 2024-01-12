/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "EmbeddedBoundary/DistanceToEB.H"
#include "WarpX.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Pusher/UpdatePosition.H"

#include <ablastr/particles/NodalFieldGather.H>

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>
#include <AMReX.H>
#include <AMReX_Algorithm.H>

struct IsOutsideDomainBoundary {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_phi;
    int m_idim;
    int m_iside;

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const SrcData& src,
                    int ip, const amrex::RandomEngine& /*engine*/) const noexcept
    {
        const auto& p = src.getSuperParticle(ip);
        if (m_iside == 0) {
            if (p.pos(m_idim) < m_plo[m_idim]) { return 1; }
        } else {
            if (p.pos(m_idim) >= m_phi[m_idim]) { return 1; }
        }
        return 0;
    }
};
#ifdef AMREX_USE_EB
struct FindBoundaryIntersection {
    const int m_index;
    const int m_step;
    const amrex::Real m_dt;
    amrex::Array4<const amrex::Real> m_phiarr;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_dxi;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator() (const DstData& dst, const SrcData& src,
                     int src_i, int dst_i) const noexcept
    {
        // Copy all particle attributes, from the source to the destination
        dst.m_aos[dst_i] = src.m_aos[src_i];
        for (int j = 0; j < SrcData::NAR; ++j) {
            dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
        }
        for (int j = 0; j < src.m_num_runtime_real; ++j) {
            dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
        }
        for (int j = 0; j < src.m_num_runtime_int; ++j) {
            dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
        }

        // Modify the position of the destination particle:
        // Move it to the point of intersection with the embedded boundary
        // (which is found by using a bisection algorithm)

        const auto& p = dst.getSuperParticle(dst_i);
        amrex::ParticleReal xp, yp, zp;
        get_particle_position( p, xp, yp, zp );
        amrex::ParticleReal const ux = dst.m_rdata[PIdx::ux][dst_i];
        amrex::ParticleReal const uy = dst.m_rdata[PIdx::uy][dst_i];
        amrex::ParticleReal const uz = dst.m_rdata[PIdx::uz][dst_i];

        // Bisection algorithm to find the point where phi(x,y,z)=0 (i.e. on the embedded boundary)

        // Temporary variables to avoid implicit capture
        amrex::Real dt = m_dt;
        amrex::Array4<const amrex::Real> phiarr = m_phiarr;
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxi = m_dxi;
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo = m_plo;

        amrex::Real dt_fraction = amrex::bisect( 0.0, 1.0,
            [=] (amrex::Real dt_frac) {
                int i, j, k;
                amrex::Real W[AMREX_SPACEDIM][2];
                amrex::Real x_temp=xp, y_temp=yp, z_temp=zp;
                UpdatePosition(x_temp, y_temp, z_temp, ux, uy, uz, -dt_frac*dt);
                ablastr::particles::compute_weights_nodal(x_temp, y_temp, z_temp, plo, dxi, i, j, k, W);
                amrex::Real phi_value  = ablastr::particles::interp_field_nodal(i, j, k, W, phiarr);
                return phi_value;
            } );

        // record the real time on the destination
        dst.m_runtime_rdata[m_index][dst_i] = m_step*m_dt + (1- dt_fraction)*m_dt;

        // Now that dt_fraction has be obtained (with bisect)
        // Save the corresponding position of the particle at the boundary
        amrex::Real x_temp=xp, y_temp=yp, z_temp=zp;
        UpdatePosition(x_temp, y_temp, z_temp, ux, uy, uz, -dt_fraction*m_dt);

        // record the components of the normal on the destination
        int i, j, k;
        amrex::Real W[AMREX_SPACEDIM][2];
        ablastr::particles::compute_weights_nodal(x_temp, y_temp, z_temp, plo, dxi, i, j, k, W);
        int ic, jc, kc; // Cell-centered indices
        amrex::Real Wc[AMREX_SPACEDIM][2]; // Cell-centered weights
        ablastr::particles::compute_weights_nodal(x_temp-0.5/dxi[0], y_temp-0.5/dxi[1], z_temp-0.5/dxi[2], plo, dxi, ic, jc, kc, Wc);
        amrex::RealVect normal = DistanceToEB::interp_normal(i, j, k, W, ic, jc, kc, Wc, phiarr, dxi);
        DistanceToEB::normalize(normal);

#if (defined WARPX_DIM_3D)
        dst.m_aos[dst_i].pos(0) = x_temp;
        dst.m_aos[dst_i].pos(1) = y_temp;
        dst.m_aos[dst_i].pos(2) = z_temp;
        //save normal components
        dst.m_runtime_rdata[m_index+1][dst_i] = normal[0];
        dst.m_runtime_rdata[m_index+2][dst_i] = normal[1];
        dst.m_runtime_rdata[m_index+3][dst_i] = normal[2];
#elif (defined WARPX_DIM_XZ)
        dst.m_aos[dst_i].pos(0) = x_temp;
        dst.m_aos[dst_i].pos(1) = z_temp;
        //save normal components
        dst.m_runtime_rdata[m_index+1][dst_i] = normal[0];
        dst.m_runtime_rdata[m_index+2][dst_i] = 0.0;
        dst.m_runtime_rdata[m_index+3][dst_i] = normal[1];
#elif (defined WARPX_DIM_RZ)
        dst.m_aos[dst_i].pos(0) = std::sqrt(x_temp*x_temp + y_temp*y_temp);
        dst.m_rdata[PIdx::theta][dst_i] = std::atan2(y_temp, x_temp);
        dst.m_aos[dst_i].pos(1) = z_temp;
        //save normal components
        amrex::Real theta=std::atan2(y_temp, x_temp);
        dst.m_runtime_rdata[m_index+1][dst_i] = normal[0]*std::cos(theta);
        dst.m_runtime_rdata[m_index+2][dst_i] = normal[0]*std::sin(theta);
        dst.m_runtime_rdata[m_index+3][dst_i] = normal[1];
#elif (defined WARPX_DIM_1D_Z)
        dst.m_aos[dst_i].pos(0) = z_temp;
        //normal not defined
        dst.m_runtime_rdata[m_index+1][dst_i] = 0.0;
        dst.m_runtime_rdata[m_index+2][dst_i] = 0.0;
        dst.m_runtime_rdata[m_index+3][dst_i] = 0.0;
#endif
    }
};
#endif

struct CopyAndTimestamp {
    int m_index;
    int m_step;
    amrex::Real m_dt;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator() (const DstData& dst, const SrcData& src,
                     int src_i, int dst_i) const noexcept
    {
        dst.m_aos[dst_i] = src.m_aos[src_i];
        for (int j = 0; j < SrcData::NAR; ++j) {
            dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
        }
        for (int j = 0; j < src.m_num_runtime_real; ++j) {
            dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
        }
        for (int j = 0; j < src.m_num_runtime_int; ++j) {
            dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
        }
        dst.m_runtime_rdata[m_index][dst_i] = m_step*m_dt;
    }
};