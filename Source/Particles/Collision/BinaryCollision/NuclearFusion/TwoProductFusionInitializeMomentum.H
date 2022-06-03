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
                            SoaData_type& soa_product1, SoaData_type& soa_product2,
                            const index_type& idx_1, const index_type& idx_2,
                            const index_type& idx_product1_start, const index_type& idx_product2_start,
                            const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
                            const amrex::ParticleReal& mp1, const amrex::ParticleReal& mp2,
                            const amrex::ParticleReal& E_fusion,
                            const amrex::RandomEngine& engine)
    {
        using namespace amrex::literals;

        constexpr amrex::ParticleReal c_sq = PhysConst::c * PhysConst::c;
        constexpr amrex::ParticleReal inv_csq = 1._prt / ( c_sq );
        // Rest energy of incident particles
        const amrex::ParticleReal E_rest_in = (m1 + m2)*c_sq;
        // Rest energy of products
        const amrex::ParticleReal E_rest_out = (mp1 + mp2)*c_sq;

        // Normalized momentum of colliding particles
        const amrex::ParticleReal u1x = soa_1.m_rdata[PIdx::ux][idx_1];
        const amrex::ParticleReal u1y = soa_1.m_rdata[PIdx::uy][idx_1];
        const amrex::ParticleReal u1z = soa_1.m_rdata[PIdx::uz][idx_1];
        const amrex::ParticleReal u2x = soa_2.m_rdata[PIdx::ux][idx_2];
        const amrex::ParticleReal u2y = soa_2.m_rdata[PIdx::uy][idx_2];
        const amrex::ParticleReal u2z = soa_2.m_rdata[PIdx::uz][idx_2];

        // Compute Lorentz factor gamma in the lab frame
        const amrex::ParticleReal g1 = std::sqrt( 1._prt + (u1x*u1x+u1y*u1y+u1z*u1z)*inv_csq );
        const amrex::ParticleReal g2 = std::sqrt( 1._prt + (u2x*u2x+u2y*u2y+u2z*u2z)*inv_csq );

        // Compute momenta
        const amrex::ParticleReal p1x = u1x * m1;
        const amrex::ParticleReal p1y = u1y * m1;
        const amrex::ParticleReal p1z = u1z * m1;
        const amrex::ParticleReal p2x = u2x * m2;
        const amrex::ParticleReal p2y = u2y * m2;
        const amrex::ParticleReal p2z = u2z * m2;
        // Square norm of the total (sum between the two particles) momenta in the lab frame
        auto constexpr pow2 = [](double const x) { return x*x; };
        const amrex::ParticleReal p_total_sq =  pow2(p1x+p2x) +
                                                pow2(p1y+p2y) +
                                                pow2(p1z+p2z);

        // Total energy of incident species in the lab frame
        const amrex::ParticleReal E_lab = (m1 * g1 + m2 * g2) * c_sq;
        // Total energy squared of proton+boron in the center of mass frame, calculated using the
        // Lorentz invariance of the four-momentum norm
        const amrex::ParticleReal E_star_sq = E_lab*E_lab - c_sq*p_total_sq;
        // Total energy squared of the products in the center of mass frame
        // In principle, the term - E_rest_in + E_rest_out + E_fusion is not needed and equal to
        // zero (i.e. the energy liberated during fusion is equal to the mass difference). However,
        // due to possible inconsistencies in how the mass is defined in the code, it is
        // probably more robust to subtract the rest masses and to add the fusion energy to the
        // total kinetic energy.
        const amrex::ParticleReal E_star_f_sq = pow2(std::sqrt(E_star_sq)
                                                         - E_rest_in + E_rest_out + E_fusion);

        // Square of the norm of the momentum of the products in the center of mass frame
        // Formula obtained by inverting E^2 = p^2*c^2 + m^2*c^4 in the COM frame for each particle
        auto constexpr pow3 = [](double const x) { return x*x*x; };
        const amrex::ParticleReal p_star_f_sq =
            E_star_f_sq*0.25_prt*inv_csq - (mp1*mp1 + mp2*mp2)*c_sq*0.5_prt +
            pow3(c_sq)*0.25_prt * pow2(mp1*mp1 - mp2*mp2) / E_star_f_sq;

        // Compute momentum of first alpha in the center of mass frame, assuming isotropic
        // distribution
        amrex::ParticleReal px_star, py_star, pz_star;
        ParticleUtils::RandomizeVelocity(px_star, py_star, pz_star, std::sqrt(p_star_f_sq),
                                         engine);

        // Next step is to convert momenta to lab frame
        amrex::ParticleReal px_out1, py_out1, pz_out1;
        // Preliminary calculation: compute center of mass velocity vc
        const amrex::ParticleReal mass_g = m1 * g1 + m2 * g2;
        const amrex::ParticleReal vcx    = (p1x+p2x) / mass_g;
        const amrex::ParticleReal vcy    = (p1y+p2y) / mass_g;
        const amrex::ParticleReal vcz    = (p1z+p2z) / mass_g;
        const amrex::ParticleReal vc_sq   = vcx*vcx + vcy*vcy + vcz*vcz;

        // Convert momentum of first alpha to lab frame, using equation (13) of F. Perez et al.,
        // Phys.Plasmas.19.083104 (2012)
        if ( vc_sq > std::numeric_limits<amrex::ParticleReal>::min() )
        {
            const amrex::ParticleReal gc = 1._prt / std::sqrt( 1._prt - vc_sq*inv_csq );
            const amrex::ParticleReal g_star = std::sqrt(1._prt + p_star_f_sq / (mp1*mp1*c_sq));
            const amrex::ParticleReal vcDps = vcx*px_star + vcy*py_star + vcz*pz_star;
            const amrex::ParticleReal factor0 = (gc-1._prt)/vc_sq;
            const amrex::ParticleReal factor = factor0*vcDps + mp1*g_star*gc;
            px_out1 = px_star + vcx * factor;
            py_out1 = py_star + vcy * factor;
            pz_out1 = pz_star + vcz * factor;
        }
        else // If center of mass velocity is zero, we are already in the lab frame
        {
            px_out1 = px_star;
            py_out1 = py_star;
            pz_out1 = pz_star;
        }

        // Compute momentum of beryllium in lab frame, using total momentum conservation
        const amrex::ParticleReal px_out2 = p1x + p2x - px_out1;
        const amrex::ParticleReal py_out2 = p1y + p2y - py_out1;
        const amrex::ParticleReal pz_out2 = p1z + p2z - pz_out1;

        // Fill momentum of product species (note that we actually
        // create 4 products, 2 at the position of each incident particle)
        soa_product1.m_rdata[PIdx::ux][idx_product1_start] = px_out1/mp1;
        soa_product1.m_rdata[PIdx::uy][idx_product1_start] = py_out1/mp1;
        soa_product1.m_rdata[PIdx::uz][idx_product1_start] = pz_out1/mp1;
        soa_product1.m_rdata[PIdx::ux][idx_product1_start + 1] = px_out1/mp1;
        soa_product1.m_rdata[PIdx::uy][idx_product1_start + 1] = py_out1/mp1;
        soa_product1.m_rdata[PIdx::uz][idx_product1_start + 1] = pz_out1/mp1;
        soa_product2.m_rdata[PIdx::ux][idx_product2_start] = px_out2/mp2;
        soa_product2.m_rdata[PIdx::uy][idx_product2_start] = py_out2/mp2;
        soa_product2.m_rdata[PIdx::uz][idx_product2_start] = pz_out2/mp2;
        soa_product2.m_rdata[PIdx::ux][idx_product2_start + 1] = px_out2/mp2;
        soa_product2.m_rdata[PIdx::uy][idx_product2_start + 1] = py_out2/mp2;
        soa_product2.m_rdata[PIdx::uz][idx_product2_start + 1] = pz_out2/mp2;
    }

}

#endif // TWO_PRODUCT_FUSION_INITIALIZE_MOMENTUM_H
