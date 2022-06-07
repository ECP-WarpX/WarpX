/* Copyright 2022 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef TWO_PRODUCT_FUSION_INITIALIZE_MOMENTUM_H
#define TWO_PRODUCT_FUSION_INITIALIZE_MOMENTUM_H

#include "TwoProductFusionUtil.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_DenseBins.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace {
    // Define shortcuts for frequently-used type names
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleBins = amrex::DenseBins<ParticleType>;
    using index_type = ParticleBins::index_type;

    /**
     * \brief This function initializes the momentum of the product particles,
     * in a fusion event where only two products are produced.
     * (In this case, conservation of energy and momentum determines
     * the amplitude of the momentum of the particles exactly.)
     * We assume that the emission of the product is isotropic in the center-of-mass frame
     *
     * @param[in] soa_1 struct of array data of the first colliding species
     * @param[in] soa_2 struct of array data of the second colliding species
     * or tritium)
     * @param[out] ...
     * @param[in] idx_1 index of first colliding macroparticle
     * @param[in] idx_2 index of second colliding macroparticle
     * @param[in]...
     * @param[in] m1 mass of first colliding species
     * @param[in] m2 mass of second colliding species
     * @param[in] engine the random engine
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void TwoProductFusionInitializeMomentum (
                            const SoaData_type& soa_1, const SoaData_type& soa_2,
                            SoaData_type& soa_product1, SoaData_type& soa_product2,
                            const index_type& idx_1, const index_type& idx_2,
                            const index_type& idx_product1_start, const index_type& idx_product2_start,
                            const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
                            const amrex::ParticleReal& mp1, const amrex::ParticleReal& mp2,
                            const amrex::ParticleReal& E_fusion,
                            const amrex::RandomEngine& engine)
    {
        using namespace amrex::literals;

        amrex::ParticleReal ux1_out = 0.0, uy1_out = 0.0, uz1_out = 0.0;
        amrex::ParticleReal ux2_out = 0.0, uy2_out = 0.0, uz2_out = 0.0;

        TwoProductFusionComputeProductMomenta(
            soa_1.m_rdata[PIdx::ux][idx_1],
            soa_1.m_rdata[PIdx::uy][idx_1],
            soa_1.m_rdata[PIdx::uz][idx_1], m1,
            soa_2.m_rdata[PIdx::ux][idx_2],
            soa_2.m_rdata[PIdx::uy][idx_2],
            soa_2.m_rdata[PIdx::uz][idx_2], m2,
            ux1_out, uy1_out, uz1_out, mp1,
            uy2_out, uy2_out, uz2_out, mp2,
            E_fusion,
            engine);

        // Fill momentum of product species (note that we actually
        // create 4 products, 2 at the position of each incident particle)
        soa_product1.m_rdata[PIdx::ux][idx_product1_start] = ux1_out;
        soa_product1.m_rdata[PIdx::uy][idx_product1_start] = uy1_out;
        soa_product1.m_rdata[PIdx::uz][idx_product1_start] = uz1_out;
        soa_product1.m_rdata[PIdx::ux][idx_product1_start + 1] = ux1_out;
        soa_product1.m_rdata[PIdx::uy][idx_product1_start + 1] = uy1_out;
        soa_product1.m_rdata[PIdx::uz][idx_product1_start + 1] = uz1_out;
        soa_product2.m_rdata[PIdx::ux][idx_product2_start] = ux2_out;
        soa_product2.m_rdata[PIdx::uy][idx_product2_start] = uy2_out;
        soa_product2.m_rdata[PIdx::uz][idx_product2_start] = uz2_out;
        soa_product2.m_rdata[PIdx::ux][idx_product2_start + 1] = ux2_out;
        soa_product2.m_rdata[PIdx::uy][idx_product2_start + 1] = uy2_out;
        soa_product2.m_rdata[PIdx::uz][idx_product2_start + 1] = uz2_out;
    }
}

#endif // TWO_PRODUCT_FUSION_INITIALIZE_MOMENTUM_H
