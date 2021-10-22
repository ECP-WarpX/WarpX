/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PROTON_BORON_FUSION_INITIALIZE_MOMENTUM_H
#define PROTON_BORON_FUSION_INITIALIZE_MOMENTUM_H

#include "Particles/WarpXParticleContainer.H"

#include <AMReX_DenseBins.H>
#include <AMReX_REAL.H>

namespace {
    // Define shortcuts for frequently-used type names
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;

    /**
     * \brief This function initializes the momentum of the alpha particles produced from
     * proton-boron fusion. It needs to be implemented and currently only zero-initializes the
     * momenta.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void ProtonBoronFusionInitializeMomentum (
                            const SoaData_type /*soa_1*/, const SoaData_type /*soa_2*/,
                            const SoaData_type soa_alpha,
                            const index_type& /*idx_1*/, const index_type& /*idx_2*/,
                            const index_type& idx_alpha_start,
                            const amrex::ParticleReal& /*m1*/, const amrex::ParticleReal& /*m2*/)
    {
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 1] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 1] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 1] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 2] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 2] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 2] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 3] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 3] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 3] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 4] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 4] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 4] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 5] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 5] = amrex::ParticleReal(0.);
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 5] = amrex::ParticleReal(0.);
    }

}

#endif // PROTON_BORON_FUSION_INITIALIZE_MOMENTUM_H
