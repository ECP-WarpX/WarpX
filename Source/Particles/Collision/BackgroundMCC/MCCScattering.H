/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_
#define WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_

#include "MCCProcess.H"

#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

/** @file
 *
 * This file contains the implementation of the scattering processes available
 * in the MCC handling.
 */

/** \brief Function to perform isotropic elastic scattering of a particle. A
 * hard sphere collision is assumed independent of the projectile velocity.
 *
 * @param[in,out] ux,uy,uz colliding particle's velocity
 * @param[in] engine random number generator.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void IsotropicScattering ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                           amrex::ParticleReal& uz, amrex::RandomEngine const& engine )
{
    using std::sqrt;

    // istropically scatter the particle
    amrex::Real const mag = sqrt(ux*ux + uy*uy + uz*uz);
    ParticleUtils::RandomizeVelocity(ux, uy, uz, mag, engine);
}


/** \brief Function to sample a scattering angle using the Wentzel-Moliere
 * model (Fernandez-Varea et al. https://doi.org/10.1016/0168-583X(93)95827-R)
 *
 * @param[in,out] cos_chi cosine of scattering angle
 * @param[in] A screening parameter (see eq. 32 of Fernandez-Varea et al.)
 * @param[in] engine random number generator.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void getWentzelScatteringAngle ( amrex::ParticleReal& cos_chi,
                                 const amrex::ParticleReal A,
                                 amrex::RandomEngine const& engine )
{
    // draw a scattering angle from the inverted cummulative
    // distribution function (see eq. 2.16 of Lim, 2007)
    cos_chi = (2.0_prt * A + 1.0_prt
        + 2.0_prt * A * (A + 1.0_prt) / (amrex::Random(engine) - A - 1.0_prt)
    );
}


/** \brief Function to perform elastic scattering of a particle in the
 * stationary frame of the target. The particle velocity is rotated so that it
 * is travelling only in the z'-direction. The Wentzel-Moliere model is used to
 * determine the scattering angle in the primed frame. The particle velocity is
 * then updated and rotated back to the unprimed frame.
 * See Fernandez-Varea et al. https://doi.org/10.1016/0168-583X(93)95827-R and
 * the PhD thesis of Chul-Hyun Lim, 2007, Berkeley, "The interaction of
 * energetic charged particles with gas and boundaries in the particle
 * simulation of plasmas".
 *
 * @param[in,out] ux,uy,uz colliding particle's velocity
 * @param[in] uCOM_x,uCOM_y,uCOM_z velocity of the center of momentum frame.
 * @param[in] engine random number generator.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ElasticScattering ( amrex::ParticleReal& vx, amrex::ParticleReal& vy,
                         amrex::ParticleReal& vz,
                         amrex::ParticleReal const E_penalty,
                         amrex::ParticleReal const m,
                         amrex::ParticleReal const M,
                         amrex::ParticleReal const A_prefactor,
                         amrex::ParticleReal const Z,
                         amrex::RandomEngine const& engine,
                         amrex::Real cos_chi=2.0_rt )
{
    // special case for back-scattering with equal masses;
    // return early to avoid unnecessary calculations
    if (M == m && cos_chi == -1.0_rt) {
        vx *= -1.0_prt;
        vy *= -1.0_prt;
        vz *= -1.0_prt;
        return;
    }

    // precalculate often used value
    constexpr auto c2 = PhysConst::c * PhysConst::c;
    const auto mc2 = m*c2;

    amrex::ParticleReal E_kin, v2, gamma;

    // calculate the remaining particle energy and from that its speed
    v2 = vx*vx + vy*vy + vz*vz;
    ParticleUtils::getEnergy(v2, m, gamma, E_kin);
    E_kin = (E_kin - E_penalty) * PhysConst::q_e;
    const auto p2 = E_kin * (E_kin + 2.0_prt*mc2) / c2;

    if (E_penalty != 0.0_prt) {
        const auto scaling_factor = sqrt(p2 / v2) / m;
        // scale the velocity vector based its new magnitude (after substracting
        // the energy penalty of the collision)
        vx *= scaling_factor;
        vy *= scaling_factor;
        vz *= scaling_factor;
        v2 = vx*vx + vy*vy + vz*vz;
    }

    // at small energies the Wentzel model gives isotropic scattering; if we are
    // in such a situation and a scattering angle was not given we can save
    // some compution
    if ((cos_chi == 2.0_prt) && (E_kin / PhysConst::q_e < 100.0_prt)) {
        IsotropicScattering(vx, vy, vz, engine);
        return;
    }

    // sample a scattering angle if one was not given
    if (cos_chi == 2.0_prt) {
        // calculate the full screening parameter (using eq. 32 of Fernandez-Varea et al.)
        auto beta = v2/c2;
        auto A = A_prefactor / p2
          * (1.13_prt + 3.76_prt * PhysConst::alpha * PhysConst::alpha
             * Z * Z / beta
        );
        getWentzelScatteringAngle(cos_chi, A, engine);
    }

    // the Wentzel model gives a scattering angle relative to the current
    // velocity direction of the particle, so we rotate the velocity vector
    // to point along z' before applying the scattering angle
    amrex::ParticleReal theta, phi, ux, uy, uz;
    int sign;
    ParticleUtils::getRotationAngles(vx, vy, vz, theta, phi, sign);
    ux = 0.0_prt;
    uy = 0.0_prt;
    uz = sign*sqrt(v2);

    // get the new velocity components in the rotated frame
    ux = uz * sqrt(1.0_prt - cos_chi*cos_chi);
    uz = uz * cos_chi;

    // choose random angle to rotate velocity around z-axis
    const auto omega = amrex::Random(engine) * 2.0_rt * MathConst::pi;
    uy = ux * sin(omega);
    ux = ux * cos(omega);

    // now rotate the velocity back to the original collision frame
    // (where the target was stationary)
    const auto sin_theta = sin(theta);
    const auto cos_theta = cos(theta);
    const auto sin_phi = sin(phi);
    const auto cos_phi = cos(phi);

    vx = ux * cos_theta * cos_phi + uy * sin_theta + uz * cos_theta * sin_phi;
    vy = -ux * sin_theta * cos_phi + uy * cos_theta - uz * sin_theta * sin_phi;
    vz = -ux * sin_phi + uz * cos_phi;
}


/**
 * \brief Filter functor for impact ionization
 */
class ImpactIonizationFilterFunc
{
public:

    /**
    * \brief Constructor of the ImpactIonizationFilterFunc functor.
    *
    * This function sample a random number and compares it to the total
    * collision probability to see if the given particle should be considered
    * for an ionization event. If the particle passes this stage the collision
    * cross-section is calculated given its energy and another random number
    * is used to determine whether it actually collides.
    *
    * @param[in] mcc_process an MCCProcess object associated with the ionization
    * @param[in] mass colliding particle's mass (could also assume electron)
    * @param[in] total_collision_prob total probability for a collision to occur
    * @param[in] nu_max maximum collision frequency
    * @param[in] n_a_func ParserExecutor<4> function to get the background
                 density in m^-3 as a function of space and time
    * @param[in] t the current simulation time
    */
    ImpactIonizationFilterFunc(
        MCCProcess const& mcc_process,
        amrex::Real const mass,
        amrex::Real const total_collision_prob,
        amrex::Real const nu_max,
        amrex::ParserExecutor<4> const& n_a_func,
        amrex::Real t
    ) : m_mcc_process(mcc_process.executor()), m_mass(mass),
        m_total_collision_prob(total_collision_prob),
        m_nu_max(nu_max), m_n_a_func(n_a_func), m_t(t) { }

    /**
    * \brief Functor call. This method determines if a given (electron) particle
    * should undergo an ionization collision.
    *
    * @param[in] ptd particle tile data
    * @param[in] i particle index
    * @param[in] engine the random number state and factory
    * @return true if a collision occurs, false otherwise
    */
    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (
        const PData& ptd, int const i, amrex::RandomEngine const& engine
    ) const noexcept
    {
        using namespace amrex;
        using std::sqrt;

        // determine if this particle should collide
        if (Random(engine) > m_total_collision_prob) return false;

        // get references to the particle to get its position
        const auto& p = ptd.getSuperParticle(i);
        ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        // calculate neutral density at particle location
        const Real n_a = m_n_a_func(x, y, z, m_t);

        // get the particle velocity
        const ParticleReal ux = ptd.m_rdata[PIdx::ux][i];
        const ParticleReal uy = ptd.m_rdata[PIdx::uy][i];
        const ParticleReal uz = ptd.m_rdata[PIdx::uz][i];

        // calculate kinetic energy
        const ParticleReal v_coll2 = ux*ux + uy*uy + uz*uz;
        const ParticleReal E_coll = 0.5_prt * m_mass * v_coll2 / PhysConst::q_e;
        const ParticleReal v_coll = sqrt(v_coll2);

        // get collision cross-section
        const Real sigma_E = m_mcc_process.getCrossSection(E_coll);

        // calculate normalized collision frequency
        const Real nu_i = n_a * sigma_E * v_coll / m_nu_max;

        // check if this collision should be performed
        return (Random(engine) <= nu_i);
    }

private:
    MCCProcess::Executor m_mcc_process;
    amrex::Real m_mass;
    amrex::Real m_total_collision_prob = 0;
    amrex::Real m_nu_max;
    amrex::ParserExecutor<4> m_n_a_func;
    amrex::Real m_t;
};


/**
 * \brief Transform functor for impact ionization
 */
class ImpactIonizationTransformFunc
{
public:

    /**
    * \brief Constructor of the ImpactIonizationTransformFunc functor.
    *
    * The transform is responsible for appropriately decreasing the kinetic
    * energy of the colliding particle and assigning appropriate velocities
    * to the two newly created particles. To this end the energy cost of
    * ionization is passed to the constructor as well as the mass of the
    * colliding species and the standard deviation of the ion velocity
    * (normalized temperature).
    *
    * @param[in] energy_cost energy cost of ionization
    * @param[in] mass1 mass of the colliding species
    * @param[in] sqrt_kb_m value of sqrt(kB/m), where kB is Boltzmann's constant
                 and m is the background neutral mass
    * @param[in] T_a_func ParserExecutor<4> function to get the background
                 temperature in Kelvin as a function of space and time
    * @param[in] t the current simulation time
    */
    ImpactIonizationTransformFunc(
        amrex::Real energy_cost, amrex::Real mass1, amrex::Real sqrt_kb_m,
        amrex::ParserExecutor<4> const& T_a_func, amrex::Real t
    ) :  m_energy_cost(energy_cost), m_mass1(mass1),
         m_sqrt_kb_m(sqrt_kb_m), m_T_a_func(T_a_func), m_t(t) { }

    /**
    * \brief Functor call. It determines the properties of the generated pair
    * and decreases the kinetic energy of the colliding particle. Inputs are
    * standard from the FilterCopyTransfrom::filterCopyTransformParticles
    * function.
    *
    * @param[in,out] dst1 target species 1 (electrons)
    * @param[in,out] dst2 target species 2 (ions)
    * @param[in] src source species (electrons)
    * @param[in] i_src particle index of the source species
    * @param[in] i_dst1 particle index of target species 1
    * @param[in] i_dst2 particle index of target species 2
    * @param[in] engine random number generator engine
    */
    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst1, DstData& dst2, SrcData& src,
        int const i_src, int const i_dst1, int const i_dst2,
        amrex::RandomEngine const& engine) const noexcept
    {
        using namespace amrex;
        using std::sqrt;

        // get references to the particle to get its position
        const auto& p = src.getSuperParticle(i_src);
        ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        // calculate standard deviation in neutral velocity distribution using
        // the local temperature
        const Real ion_vel_std = m_sqrt_kb_m * sqrt(m_T_a_func(x, y, z, m_t));

        // get references to the original particle's velocity
        auto& ux = src.m_rdata[PIdx::ux][i_src];
        auto& uy = src.m_rdata[PIdx::uy][i_src];
        auto& uz = src.m_rdata[PIdx::uz][i_src];

        // get references to the new particles' velocities
        auto& e_ux = dst1.m_rdata[PIdx::ux][i_dst1];
        auto& e_uy = dst1.m_rdata[PIdx::uy][i_dst1];
        auto& e_uz = dst1.m_rdata[PIdx::uz][i_dst1];
        auto& i_ux = dst2.m_rdata[PIdx::ux][i_dst2];
        auto& i_uy = dst2.m_rdata[PIdx::uy][i_dst2];
        auto& i_uz = dst2.m_rdata[PIdx::uz][i_dst2];

        // calculate kinetic energy
        const ParticleReal v_coll2 = ux*ux + uy*uy + uz*uz;
        const ParticleReal E_coll = 0.5_prt * m_mass1 * v_coll2 / PhysConst::q_e;

        // get the left over energy
        amrex::Real E_remaining = E_coll - m_energy_cost;

        // each electron gets half the energy (could change this later)
        amrex::Real vp = sqrt(
            2.0_prt / m_mass1 * PhysConst::q_e * E_remaining / 2.0_prt
        );

        // isotropically scatter electrons
        ParticleUtils::RandomizeVelocity(ux, uy, uz, vp, engine);
        ParticleUtils::RandomizeVelocity(e_ux, e_uy, e_uz, vp, engine);

        // get velocities for the ion from a Maxwellian distribution
        i_ux = ion_vel_std * RandomNormal(0_prt, 1.0_prt, engine);
        i_uy = ion_vel_std * RandomNormal(0_prt, 1.0_prt, engine);
        i_uz = ion_vel_std * RandomNormal(0_prt, 1.0_prt, engine);
    }

private:
    amrex::Real m_energy_cost;
    amrex::Real m_mass1;
    amrex::Real m_sqrt_kb_m;
    amrex::ParserExecutor<4> m_T_a_func;
    amrex::Real m_t;
};
#endif // WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_
