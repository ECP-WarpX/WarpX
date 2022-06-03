/* Copyright 2022 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef TWO_PRODUCT_FUSION_INITIALIZE_MOMENTUM_H
#define TWO_PRODUCT_FUSION_INITIALIZE_MOMENTUM_H

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
                            SoaData_type& soa_alpha,
                            const index_type& idx_1, const index_type& idx_2,
                            const index_type& idx_alpha_start,
                            const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
                            const amrex::RandomEngine& engine)
    {

        using namespace amrex::literals;
    }

}

#endif // TWO_PRODUCT_FUSION_INITIALIZE_MOMENTUM_H
