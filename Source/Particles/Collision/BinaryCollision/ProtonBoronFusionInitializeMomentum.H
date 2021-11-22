/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef PROTON_BORON_FUSION_INITIALIZE_MOMENTUM_H
#define PROTON_BORON_FUSION_INITIALIZE_MOMENTUM_H

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
     * \brief This function initializes the momentum of the alpha particles produced from
     * proton-boron fusion. The momentum is initialized by assuming that the fusion of a proton
     * with a boron nucleus into 3 alphas takes place in two steps. In the first step, the proton
     * and the boron fuse into a beryllium nucleus and an alpha particle. In the second step, the
     * beryllium decays into two alpha particles. The first step produces 8.59009 MeV of kinetic
     * energy while the second step produces 91.8984 keV of kinetic energy. This two-step process
     * is considered to be the dominant process of proton+boron fusion into alphas (see
     * Becker et al., Zeitschrift für Physik A Atomic Nuclei, 327(3), 341-355 (1987)).
     * For each step, we assume in this function that the particles are emitted isotropically in
     * the corresponding center of mass frame (center of mass frame of proton + boron for the
     * creation of first alpha+beryllium and rest frame of beryllium for the creation of second and
     * third alphas). This isotropic assumption is exact for the second step but is only an
     * approximation for the first step.
     *
     * @param[in] soa_1 struct of array data of the first colliding species (can be either proton
     * or boron)
     * @param[in] soa_2 struct of array data of the second colliding species (can be either proton
     * or boron)
     * @param[out] soa_alpha struct of array data of the alpha species
     * @param[in] idx_1 index of first colliding macroparticle
     * @param[in] idx_2 index of second colliding macroparticle
     * @param[in] idx_alpha_start index of first produced alpha macroparticle
     * @param[in] m1 mass of first colliding species
     * @param[in] m2 mass of second colliding species
     * @param[in] engine the random engine
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void ProtonBoronFusionInitializeMomentum (
                            const SoaData_type& soa_1, const SoaData_type& soa_2,
                            SoaData_type& soa_alpha,
                            const index_type& idx_1, const index_type& idx_2,
                            const index_type& idx_alpha_start,
                            const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
                            const amrex::RandomEngine& engine)
    {
        // General notations in this function:
        //     x_sq denotes the square of x
        //     x_star denotes the value of x in the proton+boron center of mass frame
        //     x_Bestar denotes the value of x in the Beryllium rest frame

        using namespace amrex::literals;

        constexpr amrex::ParticleReal mev_to_joule = PhysConst::q_e*1.e6_prt;
        // Energy produced in the fusion reaction proton + boron11 -> Beryllium8 + alpha
        // cf. Janis book of proton-induced cross-sections (2019)
        constexpr amrex::ParticleReal E_fusion = 8.59009_prt*mev_to_joule;
        // Energy produced when Beryllium8 decays into two alphas
        // cf. JEFF-3.3 radioactive decay data library (2017)
        constexpr amrex::ParticleReal E_decay = 0.0918984_prt*mev_to_joule;

        constexpr amrex::ParticleReal m_alpha = PhysConst::m_p * 3.97369_prt;
        constexpr amrex::ParticleReal m_beryllium = PhysConst::m_p * 7.94748_prt;
        constexpr amrex::ParticleReal c_sq = PhysConst::c * PhysConst::c;
        constexpr amrex::ParticleReal inv_csq = 1._prt / ( c_sq );
        // Rest energy of proton+boron
        const amrex::ParticleReal E_rest_pb = (m1 + m2)*c_sq;
        // Rest energy of alpha+beryllium
        constexpr amrex::ParticleReal E_rest_abe = (m_alpha + m_beryllium)*c_sq;

        // These constexprs underflow in single precision because we use SI units and sometimes the
        // code won't compile. Nuclear fusion module does not currently work with single precision
        // anyways so these #ifdef are just here to make sure that the code always compiles with
        // single precision.
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
        const amrex::ParticleReal ma_sq = m_alpha*m_alpha;
        const amrex::ParticleReal mBe_sq = m_beryllium*m_beryllium;
#else
        constexpr amrex::ParticleReal ma_sq = m_alpha*m_alpha;
        constexpr amrex::ParticleReal mBe_sq = m_beryllium*m_beryllium;
#endif

        // Normalized momentum of colliding particles (proton and boron)
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
        const amrex::ParticleReal p_total_sq =  std::pow(p1x+p2x, 2) +
                                                std::pow(p1y+p2y, 2) +
                                                std::pow(p1z+p2z, 2);

        // Total energy of proton+boron in the lab frame
        const amrex::ParticleReal E_lab = (m1 * g1 + m2 * g2) * c_sq;
        // Total energy squared of proton+boron in the center of mass frame, calculated using the
        // Lorentz invariance of the four-momentum norm
        const amrex::ParticleReal E_star_sq = E_lab*E_lab - c_sq*p_total_sq;
        // Total energy squared of beryllium+alpha in the center of mass frame
        // In principle, the term - E_rest_pb + E_rest_abe + E_fusion is not needed and equal to
        // zero (i.e. the energy liberated during fusion is equal to the mass difference). However,
        // due to possible inconsistencies in how the mass is defined in the code (e.g. currently,
        // the mass of hydrogen is the mass of the proton, not including the electron, while the
        // mass of the other elements is the atomic mass, that includes the electrons mass), it is
        // probably more robust to subtract the rest masses and to add the fusion energy to the
        // total kinetic energy.
        const amrex::ParticleReal E_star_f_sq = std::pow(std::sqrt(E_star_sq)
                                                         - E_rest_pb + E_rest_abe + E_fusion, 2);

        // Square of the norm of the momentum of Beryllium or alpha in the center of mass frame
        // Formula obtained by inverting E^2 = p^2*c^2 + m^2*c^4 in the COM frame for each particle
        const amrex::ParticleReal p_star_f_sq =
            E_star_f_sq*0.25_prt*inv_csq - (ma_sq + mBe_sq)*c_sq*0.5_prt +
            std::pow(c_sq, 3)*0.25_prt*std::pow(ma_sq - mBe_sq, 2)/E_star_f_sq;

        // Compute momentum of first alpha in the center of mass frame, assuming isotropic
        // distribution
        amrex::ParticleReal px_star, py_star, pz_star;
        ParticleUtils::RandomizeVelocity(px_star, py_star, pz_star, std::sqrt(p_star_f_sq),
                                         engine);

        // Next step is to convert momentum of first alpha to lab frame
        amrex::ParticleReal px_alpha1, py_alpha1, pz_alpha1;
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
            const amrex::ParticleReal g_star = std::sqrt(1._prt + p_star_f_sq / (ma_sq*c_sq));
            const amrex::ParticleReal vcDps = vcx*px_star + vcy*py_star + vcz*pz_star;
            const amrex::ParticleReal factor0 = (gc-1._prt)/vc_sq;
            const amrex::ParticleReal factor = factor0*vcDps + m_alpha*g_star*gc;
            px_alpha1 = px_star + vcx * factor;
            py_alpha1 = py_star + vcy * factor;
            pz_alpha1 = pz_star + vcz * factor;
        }
        else // If center of mass velocity is zero, we are already in the lab frame
        {
            px_alpha1 = px_star;
            py_alpha1 = py_star;
            pz_alpha1 = pz_star;
        }

        // Compute momentum of beryllium in lab frame, using total momentum conservation
        const amrex::ParticleReal px_Be = p1x + p2x - px_alpha1;
        const amrex::ParticleReal py_Be = p1y + p2y - py_alpha1;
        const amrex::ParticleReal pz_Be = p1z + p2z - pz_alpha1;

        // Compute momentum norm of second and third alphas in Beryllium rest frame
        // Factor 0.5 is here because each alpha only gets half of the decay energy
        constexpr amrex::ParticleReal gamma_Bestar = (1._prt + 0.5_prt*E_decay/(m_alpha*c_sq));
        constexpr amrex::ParticleReal gamma_Bestar_sq_minus_one = gamma_Bestar*gamma_Bestar - 1._prt;
        const amrex::ParticleReal p_Bestar = m_alpha*PhysConst::c*std::sqrt(gamma_Bestar_sq_minus_one);

        // Compute momentum of second alpha in Beryllium rest frame, assuming isotropic distribution
        amrex::ParticleReal px_Bestar, py_Bestar, pz_Bestar;
        ParticleUtils::RandomizeVelocity(px_Bestar, py_Bestar, pz_Bestar, p_Bestar, engine);

        // Next step is to convert momentum of second alpha to lab frame
        amrex::ParticleReal px_alpha2, py_alpha2, pz_alpha2;
        // Preliminary calculation: compute Beryllium velocity v_Be
        const amrex::ParticleReal p_Be_sq   = px_Be*px_Be + py_Be*py_Be + pz_Be*pz_Be;
        const amrex::ParticleReal g_Be      = std::sqrt(1._prt + p_Be_sq / (mBe_sq*c_sq));
        const amrex::ParticleReal mg_Be     = m_beryllium*g_Be;
        const amrex::ParticleReal v_Bex     = px_Be / mg_Be;
        const amrex::ParticleReal v_Bey     = py_Be / mg_Be;
        const amrex::ParticleReal v_Bez     = pz_Be / mg_Be;
        const amrex::ParticleReal v_Be_sq   = v_Bex*v_Bex + v_Bey*v_Bey + v_Bez*v_Bez;

        // Convert momentum of second alpha to lab frame, using equation (13) of F. Perez et al.,
        // Phys.Plasmas.19.083104 (2012)
        if ( v_Be_sq > std::numeric_limits<amrex::ParticleReal>::min() )
        {
            const amrex::ParticleReal vcDps = v_Bex*px_Bestar + v_Bey*py_Bestar + v_Bez*pz_Bestar;
            const amrex::ParticleReal factor0 = (g_Be-1._prt)/v_Be_sq;
            const amrex::ParticleReal factor = factor0*vcDps + m_alpha*gamma_Bestar*g_Be;
            px_alpha2 = px_Bestar + v_Bex * factor;
            py_alpha2 = py_Bestar + v_Bey * factor;
            pz_alpha2 = pz_Bestar + v_Bez * factor;
        }
        else // If Beryllium velocity is zero, we are already in the lab frame
        {
            px_alpha2 = px_Bestar;
            py_alpha2 = py_Bestar;
            pz_alpha2 = pz_Bestar;
        }

        // Compute momentum of third alpha in lab frame, using total momentum conservation
        const amrex::ParticleReal px_alpha3 = px_Be - px_alpha2;
        const amrex::ParticleReal py_alpha3 = py_Be - py_alpha2;
        const amrex::ParticleReal pz_alpha3 = pz_Be - pz_alpha2;

        // Fill alpha species momentum data with the computed momentum (note that we actually
        // create 6 alphas, 3 at the position of the proton and 3 at the position of the boron, so
        // each computed momentum is used twice)
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start] = px_alpha1/m_alpha;
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start] = py_alpha1/m_alpha;
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start] = pz_alpha1/m_alpha;
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 1] = px_alpha1/m_alpha;
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 1] = py_alpha1/m_alpha;
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 1] = pz_alpha1/m_alpha;
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 2] = px_alpha2/m_alpha;
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 2] = py_alpha2/m_alpha;
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 2] = pz_alpha2/m_alpha;
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 3] = px_alpha2/m_alpha;
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 3] = py_alpha2/m_alpha;
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 3] = pz_alpha2/m_alpha;
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 4] = px_alpha3/m_alpha;
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 4] = py_alpha3/m_alpha;
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 4] = pz_alpha3/m_alpha;
        soa_alpha.m_rdata[PIdx::ux][idx_alpha_start + 5] = px_alpha3/m_alpha;
        soa_alpha.m_rdata[PIdx::uy][idx_alpha_start + 5] = py_alpha3/m_alpha;
        soa_alpha.m_rdata[PIdx::uz][idx_alpha_start + 5] = pz_alpha3/m_alpha;
    }

}

#endif // PROTON_BORON_FUSION_INITIALIZE_MOMENTUM_H
