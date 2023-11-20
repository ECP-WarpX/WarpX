/* Copyright 2023 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_RECORDINGPLANEPARTICLEFUNCTOR_H_
#define WARPX_RECORDINGPLANEPARTICLEFUNCTOR_H_

#include "ComputeParticleDiagFunctor.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <AMReX.H>
#include <AMReX_AmrParticles.H>



struct RecordParticles
{
    using TmpParticles = WarpXParticleContainer::TmpParticles;

    RecordParticles ( const WarpXParIter& a_pti, TmpParticles& tmp_particle_data,
                      amrex::Real z_location, int a_offset = 0 );

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator () (const SrcData& src, int i) const noexcept
    {
        amrex::ignore_unused(src);
        amrex::ParticleReal xp, yp, zp;
        m_get_position(i, xp, yp, zp);
        int Flag = 0;
        if ( ( (zp_old[i] <= m_z_location) && (zp >= m_z_location) ) ||
             ( (zp_old[i] >= m_z_location) && (zp <= m_z_location) ))
        {
            Flag = 1;
        }
        return Flag;
    }

    GetParticlePosition m_get_position;
    amrex::Real m_z_location;
    amrex::ParticleReal* AMREX_RESTRICT zp_old = nullptr;
};

class
RecordingPlaneParticleFunctor final : public ComputeParticleDiagFunctor
{
public:
    RecordingPlaneParticleFunctor (WarpXParticleContainer *pc_src, std::string species_name, int num_station_buffers);

    void operator () (PinnedMemoryParticleContainer& pc_dst, int &TotalParticleCounter, int i_buffer) const override;

    void PrepareFunctorData (const int i_buffer, bool record_particles, amrex::Real z_location, amrex::Real current_z_boost = 0., amrex::Real tlab = 0., int snapshot_full = 0) override;

private:

    WarpXParticleContainer* m_pc_src = nullptr;
    std::string m_species_name;
    int m_num_station_buffers;
    amrex::Real m_z_location;
    bool m_record_particles;
};

struct PlaneCrossingTime
{
    using TmpParticles = WarpXParticleContainer::TmpParticles;

    PlaneCrossingTime (const WarpXParIter& a_pti, TmpParticles& tmp_particle_data,
                       amrex::Real current_time,
                       amrex::Real z_station_location, int a_index, int a_offset = 0);

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator () (const DstData& dst, const SrcData& src, int i_src, int i_dst) const noexcept
    {

        amrex::ignore_unused(src);
        using namespace amrex::literals;
        amrex::ParticleReal xpnew, ypnew, zpnew;
        m_get_position(i_src, xpnew, ypnew, zpnew);

        // Copy runtime data
        for (int j = 0; j < src.m_num_runtime_real; ++j) {
            dst.m_runtime_rdata[j][i_dst] = src.m_runtime_rdata[j][i_src];
        }
        for (int j = 0; j < src.m_num_runtime_int; ++j) {
            dst.m_runtime_idata[j][i_dst] = src.m_runtime_idata[j][i_src];
        }

        // Plane crossing time
        amrex::Real gamma = std::sqrt(1._rt + m_inv_c2 *
                            ( m_uxnew[i_src] * m_uxnew[i_src]
                            + m_uynew[i_src] * m_uynew[i_src]
                            + m_uznew[i_src] * m_uznew[i_src]));
        amrex::Real vz = m_uznew[i_src] / gamma;
        amrex::Real time_to_plane = amrex::Math::abs(zpnew - m_z_station)/vz;
        dst.m_runtime_rdata[m_index][i_dst] = m_current_time - time_to_plane;

        // Interpolate attribs to get values at recording plane
        const amrex::Real weight_old = (zpnew - m_z_station)
                                     / (zpnew - m_zpold[i_src]);
        const amrex::Real weight_new = (m_z_station - m_zpold[i_src])
                                     / (zpnew - m_zpold[i_src]);
        const amrex::ParticleReal xp = m_xpold[i_src] * weight_old + xpnew * weight_new;
        const amrex::ParticleReal yp = m_ypold[i_src] * weight_old + ypnew * weight_new;
        const amrex::ParticleReal zp = m_z_station;
        const amrex::ParticleReal uxp = m_uxpold[i_src] * weight_old + m_uxnew[i_src] * weight_new;
        const amrex::ParticleReal uyp = m_uypold[i_src] * weight_old + m_uynew[i_src] * weight_new;
        const amrex::ParticleReal uzp = m_uzpold[i_src] * weight_old + m_uznew[i_src] * weight_new;
#if defined (WARPX_DIM_3D)
        dst.m_aos[i_dst].pos(0) = xp;
        dst.m_aos[i_dst].pos(1) = yp;
        dst.m_aos[i_dst].pos(2) = zp;
#elif defined (WARPX_DIM_RZ)
        dst.m_aos[i_dst].pos(0) = std::sqrt(xp*xp + yp*yp);
        dst.m_aos[i_dst].pos(1) = zp;
        dst.m_rdata[PIdx::theta][i_dst] = std::atan2(yp, xp);
#elif defined (WARPX_DIM_XZ)
        dst.m_aos[i_dst].pos(0) = xp;
        dst.m_aos[i_dst].pos(1) = zp;
        amrex::ignore_unused(yp);
#elif defined (WARPX_DIM_1D)
        dst.m_aos[i_dst].pos(0) = zp;
        amrex::ignore_unused(xp, yp);
#else
        amrex::ignore_unused(xp, yp, zp);
#endif
        dst.m_rdata[PIdx::w][i_dst] = m_wpnew[i_src];
        dst.m_rdata[PIdx::ux][i_dst] = uxp;
        dst.m_rdata[PIdx::uy][i_dst] = uyp;
        dst.m_rdata[PIdx::uz][i_dst] = uzp;

    }

    GetParticlePosition m_get_position;

    amrex::ParticleReal* AMREX_RESTRICT m_xpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_ypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_zpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_uxpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_uypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_uzpold = nullptr;

    const amrex::ParticleReal* AMREX_RESTRICT m_uxnew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uynew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uznew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_wpnew = nullptr;
    amrex::Real m_z_station;
    amrex::Real m_current_time;
    int m_index;
    amrex::Real m_Phys_c;
    amrex::Real m_inv_c2;

};
#endif
