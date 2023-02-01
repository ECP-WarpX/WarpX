/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef SINGLE_PHOTON_PHOTON_COLLISION_EVENT_H_
#define SINGLE_PHOTON_PHOTON_COLLISION_EVENT_H_

#include "LinearBreitWheelerCrossSection.H"

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
 * @param[in] m1,m2 masses NO
 * @param[in] w1,w2 effective weight of the colliding particles
 * @param[in] dt is the time step length between two collision calls.
 * @param[in] dV is the volume of the corresponding cell.
 * @param[in] pair_index is the index of the colliding pair
 * @param[out] p_mask is a mask that will be set to true if fusion occurs for that pair
 * @param[out] p_pair_reaction_weight stores the weight of the product particles
 * @param[in] event_multiplier factor used to increase the number of fusion events by
 * decreasing the weight of the produced particles
 * @param[in] multiplier_ratio factor used to take into account unsampled pairs (i.e. the fact
 * that a particle only collides with one or few particles of the other species)
 * @param[in] probability_threshold probability threshold above which we decrease the fusion
 * multiplier
 * @param[in] probability_target_value if the probability threshold is exceeded, this is used
 * to determine by how much the fusion multiplier is reduced
 * @param[in] fusion_type the physical fusion process to model NO
 * @param[in] engine the random engine.
 */
template <typename index_type>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void SinglePhotonPhotonCollisionEvent (const amrex::ParticleReal& u1x, const amrex::ParticleReal& u1y,
                               const amrex::ParticleReal& u1z, const amrex::ParticleReal& u2x,
                               const amrex::ParticleReal& u2y, const amrex::ParticleReal& u2z,
                               amrex::ParticleReal w1, amrex::ParticleReal w2,
                               const amrex::Real& dt, const amrex::ParticleReal& dV, const int& pair_index,
                               index_type* AMREX_RESTRICT p_mask,
                               amrex::ParticleReal* AMREX_RESTRICT p_pair_reaction_weight,
                               const amrex::ParticleReal& event_multiplier,
                               const int& multiplier_ratio,
                               const amrex::ParticleReal& probability_threshold,
                               const amrex::ParticleReal& probability_target_value,
                               const amrex::RandomEngine& engine)
{
    // General notations in this function:
    //     x_sq denotes the square of x
    //     x_star denotes the value of x in the center of mass frame

    using namespace amrex::literals;

    const amrex::ParticleReal w_min = amrex::min(w1, w2);
    const amrex::ParticleReal w_max = amrex::max(w1, w2);

    constexpr auto one_pr = amrex::ParticleReal(1.);
    constexpr auto two_pr = amrex::ParticleReal(2.);
    constexpr auto one_half_pr = amrex::ParticleReal(1./2.);
    constexpr amrex::ParticleReal c = PhysConst::c;
    constexpr amrex::ParticleReal c_sq = c*c;
    constexpr amrex::ParticleReal me = PhysConst::m_e;
    auto constexpr pow2 = [](double const x) { return x*x; };
  
    // Compute momenta in the lab frame
    const amrex::ParticleReal p1x = u1x * me;
    const amrex::ParticleReal p1y = u1y * me;
    const amrex::ParticleReal p1z = u1z * me;
    const amrex::ParticleReal p2x = u2x * me;
    const amrex::ParticleReal p2y = u2y * me;
    const amrex::ParticleReal p2z = u2z * me;
    const amrex::ParticleReal p1 = std::sqrt(pow2(p1x)+pow2(p1y)+pow2(p1z));
    const amrex::ParticleReal p2 = std::sqrt(pow2(p2x)+pow2(p2y)+pow2(p2z));

    // Compute Lorentz factor gamma in the lab frame 
    const amrex::ParticleReal g1 = p1 / ( me * c );
    const amrex::ParticleReal g2 = p2 / ( me * c );

    // Compute angle between photon momenta in the lab frame
    const amrex::ParticleReal cos_ang = (p1x*p2x+p1y*p2y+p1z*p2z)/(p1*p2);
    
    // Energy squared of each of the particles (both photons, electron and positron) 
    // in the center of mass frame (they are all equal),
    // calculated using the Lorentz invariance of the total four-momentum norm
    const amrex::ParticleReal E_star_sq = one_half_pr*c_sq*p1*p2*(one_pr - cos_ang);

    //amrex::Print() << "energy with gamma = " << g1 * me * c_sq << ", ...energy with p = " << p1 * PhysConst::c << "\n";
    //amrex::Print() << "energy with gamma = " << g2 * me * c_sq << ", ...energy with p = " << p2 * PhysConst::c << "\n";

    // Kinetic energy of each particle in the center of mass frame 
    const amrex::ParticleReal E_star = std::sqrt(E_star_sq);

    // Compute linear Breit Wheeler cross section as a function of kinetic energy in the center of mass frame
    auto lbw_cross_section = amrex::ParticleReal(0.);

    lbw_cross_section = LinearBreitWheelerCrossSection(E_star);

    // Lorentz factors in the center of mass frame
    const amrex::ParticleReal g1_star = E_star / (me*c_sq); 
    const amrex::ParticleReal g2_star = g1_star; 

    // relative velocity in the center of mass frame
    const amrex::ParticleReal v_rel = two_pr*c; 

    // Fusion cross section and relative velocity are computed in the center of mass frame.
    // On the other hand, the particle densities (weight over volume) in the lab frame are used. To
    // take into account this discrepancy, we need to multiply the fusion probability by the ratio
    // between the Lorentz factors in the COM frame and the Lorentz factors in the lab frame
    // (see Perez et al., Phys.Plasmas.19.083104 (2012))
    const amrex::ParticleReal lab_to_COM_factor = g1_star*g2_star/(g1*g2);

    // First estimate of probability to have fusion reaction
    amrex::ParticleReal probability_estimate = multiplier_ratio * event_multiplier *
                                lab_to_COM_factor * w_max * lbw_cross_section * v_rel * dt / dV;

    // Effective fusion multiplier
    amrex::ParticleReal event_multiplier_eff = event_multiplier;

    //amrex::Print() << " multiplier = " << event_multiplier << "\n";

    //amrex::Print() << " lab to com factor = " << lab_to_COM_factor << "\n";
    //amrex::Print() << "probability_threshold = " << probability_threshold << "::: probability_estimate = " << probability_estimate << "\n";


    // If the fusion probability is too high and the fusion multiplier greater than one, we risk to
    // systematically underestimate the fusion yield. In this case, we reduce the fusion multiplier
    // to reduce the fusion probability
    if (probability_estimate > probability_threshold)
    {
        // We aim for a fusion probability of probability_target_value but take into account
        // the constraint that the event_multiplier cannot be smaller than one
        event_multiplier_eff  = amrex::max(event_multiplier *
                                         probability_target_value / probability_estimate , one_pr);
        probability_estimate *= event_multiplier_eff/event_multiplier;
    }

    //amrex::Print() << "probability_estimate = " << probability_estimate << "\n";

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
    //amrex::Print() << "cross section = " << lbw_cross_section << ".... randon number = " << random_number << "::: probability = " << probability << "\n";

    // If we have a fusion event, set the mask the true and fill the product weight array
    if (random_number < probability)
    {
        p_mask[pair_index] = true;
        p_pair_reaction_weight[pair_index] = w_min/event_multiplier_eff;
        //amrex::Print() << "product weight = " << p_pair_reaction_weight[pair_index] << "\n";


    }
    else
    {
        p_mask[pair_index] = false;
    }

}


#endif // SINGLE_PHOTON_PHOTON_COLLISION_EVENT_H_
