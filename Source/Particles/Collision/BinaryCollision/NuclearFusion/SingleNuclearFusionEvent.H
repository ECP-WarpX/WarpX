/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef SINGLE_NUCLEAR_FUSION_EVENT_H_
#define SINGLE_NUCLEAR_FUSION_EVENT_H_

#include "ProtonBoronFusionCrossSection.H"

#include "Particles/Collision/BinaryCollision/BinaryCollisionUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>


/**
 * \brief This function computes whether the collision between two particles result in a
 * nuclear fusion event, using the algorithm described in Higginson et al., Journal of
 * Computational Physics 388, 439-453 (2019). If nuclear fusion occurs, the mask is set to true
 * for that given pair of particles and the weight of the produced particles is stored in
 * p_pair_reaction_weight.
 *
 * @tparam index_type type of the index argument
 * @param[in] u1x,u1y,u1z momenta of the first colliding particle
 * @param[in] u2x,u2y,u2z momenta of the second colliding particle
 * @param[in] m1,m2 masses
 * @param[in] w1,w2 effective weight of the colliding particles
 * @param[in] dt is the time step length between two collision calls.
 * @param[in] dV is the volume of the corresponding cell.
 * @param[in] pair_index is the index of the colliding pair
 * @param[out] p_mask is a mask that will be set to true if fusion occurs for that pair
 * @param[out] p_pair_reaction_weight stores the weight of the product particles
 * @param[in] fusion_multiplier factor used to increase the number of fusion events by
 * decreasing the weight of the produced particles
 * @param[in] multiplier_ratio factor used to take into account unsampled pairs (i.e. the fact
 * that a particle only collides with one or few particles of the other species)
 * @param[in] probability_threshold probability threshold above which we decrease the fusion
 * multiplier
 * @param[in] probability_target_value if the probability threshold is exceeded, this is used
 * to determine by how much the fusion multiplier is reduced
 * @param[in] fusion_type the physical fusion process to model
 * @param[in] engine the random engine.
 */
template <typename index_type>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void SingleNuclearFusionEvent (const amrex::ParticleReal& u1x, const amrex::ParticleReal& u1y,
                               const amrex::ParticleReal& u1z, const amrex::ParticleReal& u2x,
                               const amrex::ParticleReal& u2y, const amrex::ParticleReal& u2z,
                               const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
                               amrex::ParticleReal w1, amrex::ParticleReal w2,
                               const amrex::Real& dt, const amrex::Real& dV, const int& pair_index,
                               index_type* AMREX_RESTRICT p_mask,
                               amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight,
                               const amrex::Real& fusion_multiplier,
                               const int& multiplier_ratio,
                               const amrex::Real& probability_threshold,
                               const amrex::Real& probability_target_value,
                               const NuclearFusionType& fusion_type,
                               const amrex::RandomEngine& engine)
{
    // General notations in this function:
    //     x_sq denotes the square of x
    //     x_star denotes the value of x in the center of mass frame

    const amrex::ParticleReal w_min = amrex::min(w1, w2);
    const amrex::ParticleReal w_max = amrex::max(w1, w2);

    constexpr auto one_pr = amrex::ParticleReal(1.);
    constexpr auto inv_two_pr = amrex::ParticleReal(1./2.);
    constexpr auto inv_four_pr = amrex::ParticleReal(1./4.);
    constexpr amrex::ParticleReal c_sq = PhysConst::c * PhysConst::c;
    constexpr amrex::ParticleReal inv_csq = one_pr / ( c_sq );

    const amrex::ParticleReal m1_sq = m1*m1;
    const amrex::ParticleReal m2_sq = m2*m2;

    // Compute Lorentz factor gamma in the lab frame
    const amrex::ParticleReal g1 = std::sqrt( one_pr + (u1x*u1x+u1y*u1y+u1z*u1z)*inv_csq );
    const amrex::ParticleReal g2 = std::sqrt( one_pr + (u2x*u2x+u2y*u2y+u2z*u2z)*inv_csq );

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

    // Total energy in the lab frame
    const amrex::ParticleReal E_lab = (m1 * g1 + m2 * g2) * c_sq;
    // Total energy squared in the center of mass frame, calculated using the Lorentz invariance
    // of the four-momentum norm
    const amrex::ParticleReal E_star_sq = E_lab*E_lab - c_sq*p_total_sq;

    // Kinetic energy in the center of mass frame
    const amrex::ParticleReal E_kin_star = std::sqrt(E_star_sq) - (m1 + m2)*c_sq;

    // Compute fusion cross section as a function of kinetic energy in the center of mass frame
    auto fusion_cross_section = amrex::ParticleReal(0.);
    if (fusion_type == NuclearFusionType::ProtonBoron)
    {
        fusion_cross_section = ProtonBoronFusionCrossSection(E_kin_star);
    }

    // Square of the norm of the momentum of one of the particles in the center of mass frame
    // Formula obtained by inverting E^2 = p^2*c^2 + m^2*c^4 in the COM frame for each particle
    const amrex::ParticleReal p_star_sq =
            E_star_sq*inv_four_pr*inv_csq - (m1_sq + m2_sq)*c_sq*inv_two_pr +
            std::pow(c_sq, 3)*inv_four_pr*std::pow(m1_sq - m2_sq, 2)/E_star_sq;

    // Lorentz factors in the center of mass frame
    const amrex::ParticleReal g1_star = std::sqrt(1 + p_star_sq / (m1_sq*c_sq));
    const amrex::ParticleReal g2_star = std::sqrt(1 + p_star_sq / (m2_sq*c_sq));

    // relative velocity in the center of mass frame
    const amrex::ParticleReal v_rel = std::sqrt(p_star_sq) * (one_pr/(m1*g1_star) +
                                                              one_pr/(m2*g2_star));

    // Fusion cross section and relative velocity are computed in the center of mass frame.
    // On the other hand, the particle densities (weight over volume) in the lab frame are used. To
    // take into account this discrepancy, we need to multiply the fusion probability by the ratio
    // between the Lorentz factors in the COM frame and the Lorentz factors in the lab frame
    // (see Perez et al., Phys.Plasmas.19.083104 (2012))
    const amrex::ParticleReal lab_to_COM_factor = g1_star*g2_star/(g1*g2);

    // First estimate of probability to have fusion reaction
    amrex::ParticleReal probability_estimate = multiplier_ratio * fusion_multiplier *
                                lab_to_COM_factor * w_max * fusion_cross_section * v_rel * dt / dV;

    // Effective fusion multiplier
    amrex::ParticleReal fusion_multiplier_eff = fusion_multiplier;

    // If the fusion probability is too high and the fusion multiplier greater than one, we risk to
    // systematically underestimate the fusion yield. In this case, we reduce the fusion multiplier
    // to reduce the fusion probability
    if (probability_estimate > probability_threshold)
    {
        // We aim for a fusion probability of probability_target_value but take into account
        // the constraint that the fusion_multiplier cannot be smaller than one
        fusion_multiplier_eff  = amrex::max(fusion_multiplier *
                                         probability_target_value / probability_estimate , one_pr);
        probability_estimate *= fusion_multiplier_eff/fusion_multiplier;
    }

    // Compute actual fusion probability that is always between zero and one
    // In principle this is obtained by computing 1 - exp(-probability_estimate)
    // However, the computation of this quantity can fail numerically when probability_estimate is
    // too small (e.g. exp(-probability_estimate) returns 1 and the computation returns 0).
    // In this case, we simply use "probability_estimate" instead of 1 - exp(-probability_estimate)
    // The threshold exp_threshold at which we switch between the two formulas is determined by the
    // fact that computing the exponential is only useful if it can resolve the x^2/2 term of its
    // Taylor expansion, i.e. the square of probability_estimate should be greater than the
    // machine epsilon.
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
    constexpr auto exp_threshold = amrex::ParticleReal(1.e-3);
#else
    constexpr auto exp_threshold = amrex::ParticleReal(5.e-8);
#endif
    const amrex::ParticleReal probability = (probability_estimate < exp_threshold) ?
                                    probability_estimate: one_pr - std::exp(-probability_estimate);

    // Get a random number
    amrex::ParticleReal random_number = amrex::Random(engine);

    // If we have a fusion event, set the mask the true and fill the product weight array
    if (random_number < probability)
    {
        p_mask[pair_index] = true;
        p_pair_reaction_weight[pair_index] = w_min/fusion_multiplier_eff;
    }
    else
    {
        p_mask[pair_index] = false;
    }

}


#endif // SINGLE_NUCLEAR_FUSION_EVENT_H_
