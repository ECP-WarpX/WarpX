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
#endif

