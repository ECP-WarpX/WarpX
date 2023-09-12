/* Copyright 2023 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_STATIONPARTICLEFUNCTOR_H_
#define WARPX_STATIONPARTICLEFUNCTOR_H_

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
StationParticleFunctor final : public ComputeParticleDiagFunctor
{
public:
    StationParticleFunctor (WarpXParticleContainer *pc_src, std::string species_name, int num_station_buffers);

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
    PlaneCrossingTime (const WarpXParIter& a_pti, amrex::Real current_time,
                       amrex::Real z_station_location, int m_index, int a_offset = 0);

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator () (const DstData& dst, const SrcData& src, int i_src, int i_dst) const noexcept
    {

        amrex::ignore_unused(src, i_src);
        using namespace amrex::literals;
        amrex::ParticleReal xpnew, ypnew, zpnew;
        m_get_position(i_src, xpnew, ypnew, zpnew);

        amrex::Real gamma = std::sqrt(1._rt + m_inv_c2 *
                            ( m_uxnew[i_src] * m_uxnew[i_src]
                            + m_uynew[i_src] * m_uynew[i_src]
                            + m_uznew[i_src] * m_uznew[i_src]));

        amrex::Real vz = m_uznew[i_src] * m_Phys_c / gamma;

        amrex::Real time_to_plane = amrex::Math::abs(zpnew - m_z_station)/vz;
        dst.m_runtime_rdata[m_index][i_dst] = m_current_time - time_to_plane;
        // z-position of particles that cross the plane are set equal to the location of the plane
        dst.m_aos[i_dst].pos(2) = m_z_station;
    }

    GetParticlePosition m_get_position;
    const amrex::ParticleReal* AMREX_RESTRICT m_uxnew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uynew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uznew = nullptr;
    amrex::Real m_z_station;
    amrex::Real m_current_time;
    int m_index;
    amrex::Real m_Phys_c;
    amrex::Real m_inv_c2;

};
#endif
