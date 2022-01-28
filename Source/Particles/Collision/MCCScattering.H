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

/** \brief Function to perform elastic scattering of a particle in the lab
 * frame. The particle velocities transformed to the COM frame where a hard
 * sphere collision occurs. The resulting particle velocities are transformed
 * back to the lab frame and the input particle's velocity is updated.
 *
 * @param[in,out] ux,uy,uz colliding particle's velocity
 * @param[in] uCOM_x,uCOM_y,uCOM_z velocity of the center of momentum frame.
 * @param[in] engine random number generator.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ElasticScattering ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                         amrex::ParticleReal& uz, amrex::ParticleReal uCOM_x,
                         amrex::ParticleReal uCOM_y, amrex::ParticleReal uCOM_z,
                         amrex::RandomEngine const& engine )
{
    using std::sqrt;

    // transform to center of momentum frame
    ux -= uCOM_x;
    uy -= uCOM_y;
    uz -= uCOM_z;

    // istropically scatter the particle
    amrex::Real const mag = sqrt(ux*ux + uy*uy + uz*uz);
    ParticleUtils::RandomizeVelocity(ux, uy, uz, mag, engine);

    // transform back to lab frame
    ux += uCOM_x;
    uy += uCOM_y;
    uz += uCOM_z;
}


/** Function to perform back scattering of a particle in the lab frame.
 *
 * The particle velocity is transformed to the COM frame where it is
 * reversed. The resulting particle velocities are then transformed back to the
 * lab frame and the input particle's velocity is updated.
 *
 * @param[in,out] ux,uy,uz colliding particle's velocity
 * @param[in] uCOM_x,uCOM_y,uCOM_z velocity of the center of momentum frame.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void BackScattering ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                      amrex::ParticleReal& uz,
                      const amrex::ParticleReal uCOM_x,
                      const amrex::ParticleReal uCOM_y,
                      const amrex::ParticleReal uCOM_z)
{
    // transform to COM frame, reverse particle velocity and transform back
    using namespace amrex::literals;
    ux = -1.0_prt * ux + 2.0_prt * uCOM_x;
    uy = -1.0_prt * uy + 2.0_prt * uCOM_y;
    uz = -1.0_prt * uz + 2.0_prt * uCOM_z;
}


/** Function to perform charge exchange of an ion with a neutral particle.
 *
 * @param[in,out] ux,uy,uz colliding particle's velocity
 * @param[in] ua_x,ua_y,ua_z velocity of the neutral particle.
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ChargeExchange ( amrex::ParticleReal& ux, amrex::ParticleReal& uy,
                      amrex::ParticleReal& uz,
                      const amrex::ParticleReal ua_x,
                      const amrex::ParticleReal ua_y,
                      const amrex::ParticleReal ua_z)
{
    // swap ion velocity for neutral velocity
    ux = ua_x;
    uy = ua_y;
    uz = ua_z;
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
    * @param[in] nu_max_norm normalized maximum collision frequency, normalized
    *            by the neutral density.
    */
    ImpactIonizationFilterFunc(
        MCCProcess const& mcc_process,
        amrex::Real const mass,
        amrex::Real const total_collision_prob,
        amrex::Real const nu_max_norm
    ) : m_mcc_process(mcc_process.executor()), m_mass(mass),
        m_total_collision_prob(total_collision_prob),
        m_nu_max_norm(nu_max_norm){ }

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
        const Real nu_i = sigma_E * v_coll / m_nu_max_norm;

        // check if this collision should be performed
        return (Random(engine) <= nu_i);
    }

private:
    MCCProcess::Executor m_mcc_process;
    amrex::Real m_mass;
    amrex::Real m_total_collision_prob = 0;
    amrex::Real m_nu_max_norm;
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
    * @param[in] standard deviation of ion velocity distribution for each velocity
                 direction; for a Maxwellian distribution it should equal sqrt(kB*T/m)
    */
    ImpactIonizationTransformFunc(
        amrex::Real energy_cost, amrex::Real mass1, amrex::Real ion_vel_std
    ) :  m_energy_cost(energy_cost), m_mass1(mass1),
         m_ion_vel_std(ion_vel_std) { }

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
        i_ux = m_ion_vel_std * RandomNormal(0_prt, 1.0_prt, engine);
        i_uy = m_ion_vel_std * RandomNormal(0_prt, 1.0_prt, engine);
        i_uz = m_ion_vel_std * RandomNormal(0_prt, 1.0_prt, engine);
    }

private:
    amrex::Real m_energy_cost;
    amrex::Real m_mass1;
    amrex::Real m_ion_vel_std;
};
#endif // WARPX_PARTICLES_COLLISION_MCC_SCATTERING_H_
