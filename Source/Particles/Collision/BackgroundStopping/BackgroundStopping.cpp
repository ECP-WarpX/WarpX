/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundStopping.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

BackgroundStopping::BackgroundStopping (std::string const collision_name)
    : CollisionBase(collision_name)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() == 1,
                                     "Background stopping must have exactly one species.");

    amrex::ParmParse pp_collision_name(collision_name);

    getWithParser(pp_collision_name, "background_density", m_background_density);
    getWithParser(pp_collision_name, "background_temperature", m_background_temperature);
    getWithParser(pp_collision_name, "background_mass", m_background_mass);

}

void
BackgroundStopping::doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc)
{
    WARPX_PROFILE("BackgroundStopping::doCollisions()");
    using namespace amrex::literals;

    const amrex::Real dt = WarpX::GetInstance().getdt(0);
    if (int(std::floor(cur_time/dt)) % m_ndt != 0) return;

    auto& species = mypc->GetParticleContainerFromName(m_species_names[0]);
    amrex::Real species_mass = species.getMass();
    amrex::Real species_charge = species.getCharge();

    // Loop over refinement levels
    auto const flvl = species.finestLevel();
    for (int lev = 0; lev <= flvl; ++lev) {

        auto cost = WarpX::getCosts(lev);

        // loop over particles box by box
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(species, lev); pti.isValid(); ++pti) {
            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            amrex::Real wt = amrex::second();

            doBackgroundStoppingWithinTile(pti, dt, species_mass, species_charge);

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
            }
        }

    }
}

void BackgroundStopping::doBackgroundStoppingWithinTile (WarpXParIter& pti, amrex::Real dt,
                                                         amrex::Real species_mass, amrex::Real species_charge)
{
    using namespace amrex::literals;

    // So that CUDA code gets its intrinsic, not the host-only C++ library version
    using std::sqrt, std::abs, std::log, std::exp;

    // get particle count
    const long np = pti.numParticles();

    // get neutral properties
    amrex::Real n_a = m_background_density;
    amrex::Real T_a = m_background_temperature*PhysConst::q_e; // Converted from eV to Joules
    amrex::Real mass_a = m_background_mass;

    // get Struct-Of-Array particle data, also called attribs
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    amrex::ParallelFor(np,
        [=] AMREX_GPU_HOST_DEVICE (long ip)
        {

            // This implements the equation 14.12 from Introduction to Plasma Physics,
            // Goldston and Rutherford, the slowing down of beam ions due to collisions with electrons.
            // The equation is written as dV/dt = -alpha*V, and integrated to
            // give V(t+dt) = V(t)*exp(-alpha*dt)

            amrex::Real constexpr pi = MathConst::pi;
            amrex::Real constexpr ep0 = PhysConst::ep0;
            amrex::Real constexpr q_e = PhysConst::q_e;
            amrex::Real constexpr q_e2 = q_e*q_e;
            amrex::Real constexpr ep02 = ep0*ep0;

            amrex::Real const Zb = abs(species_charge/q_e);

            amrex::Real const vth = sqrt(3._rt*T_a/mass_a);
            amrex::Real const wp = sqrt(n_a*q_e2/(ep0*mass_a));
            amrex::Real const lambdadb = vth/wp;
            amrex::Real const lambdadb3 = lambdadb*lambdadb*lambdadb;
            amrex::Real const loglambda = log((12._rt*pi/Zb)*(n_a*lambdadb3));

            amrex::Real const pi32 = pi*sqrt(pi);
            amrex::Real const q2 = species_charge*species_charge;
            amrex::Real const T32 = T_a*sqrt(T_a);

            amrex::Real alpha = sqrt(2._rt)*n_a*q2*q_e2*sqrt(mass_a)*loglambda/(12._rt*pi32*ep02*species_mass*T32);

            ux[ip] *= exp(-alpha*dt);
            uy[ip] *= exp(-alpha*dt);
            uz[ip] *= exp(-alpha*dt);

        }
        );
}

