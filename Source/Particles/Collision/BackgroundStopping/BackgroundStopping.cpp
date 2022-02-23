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

    std::string background_density_str;
    std::string background_temperature_str;

    if (queryWithParser(pp_collision_name, "background_density", m_background_density)) {
        m_background_density_is_parsed = false;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_background_density > 0,
                 "For background stopping, the background density must be greater than 0");
    } else if (pp_collision_name.query("background_density(x,y,z,t)", background_density_str)) {
        m_background_density_is_parsed = true;
        m_background_density_parser = makeParser(background_density_str, {"x", "y", "z", "t"});
        m_background_density_func = m_background_density_parser.compile<4>();
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                 "For background stopping, the background density must be specified.");
    }

    if (queryWithParser(pp_collision_name, "background_temperature", m_background_temperature)) {
        m_background_temperature_is_parsed = false;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_background_temperature > 0,
                 "For background stopping, the background temperature must be greater than 0");
    } else if (pp_collision_name.query("background_temperature(x,y,z,t)", background_temperature_str)) {
        m_background_temperature_is_parsed = true;
        m_background_temperature_parser = makeParser(background_temperature_str, {"x", "y", "z", "t"});
        m_background_temperature_func = m_background_temperature_parser.compile<4>();
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                 "For background stopping, the background temperature must be specified.");
    }

    m_background_mass = PhysConst::m_e;
    if (queryWithParser(pp_collision_name, "background_mass", m_background_mass)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_background_mass > 0,
                 "For background stopping, the background mass must be greater than 0");
    }

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

            doBackgroundStoppingWithinTile(pti, dt, cur_time, species_mass, species_charge);

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
            }
        }

    }
}

void BackgroundStopping::doBackgroundStoppingWithinTile (WarpXParIter& pti, amrex::Real dt, amrex::Real t,
                                                         amrex::Real species_mass, amrex::Real species_charge)
{
    using namespace amrex::literals;

    // So that CUDA code gets its intrinsic, not the host-only C++ library version
    using std::sqrt, std::abs, std::log, std::exp;

    // get particle count
    const long np = pti.numParticles();

    // get background properties
    amrex::Real n_e_constant = m_background_density;
    amrex::Real T_e_constant = m_background_temperature*PhysConst::kb; // Converted from K to Joules
    amrex::Real mass_e = m_background_mass;

    // setup parsers for the background density and temperature
    bool background_density_is_parsed = m_background_density_is_parsed;
    bool background_temperature_is_parsed = m_background_temperature_is_parsed;
    auto n_e_func = m_background_density_func;
    auto T_e_func = m_background_temperature_func;

    // get Struct-Of-Array particle data, also called attribs
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    // May be needed to evaluate the density and/or temperature functions
    auto GetPosition = GetParticlePosition(pti);

    amrex::ParallelFor(np,
        [=] AMREX_GPU_HOST_DEVICE (long ip)
        {

            amrex::Real n_e = n_e_constant;
            amrex::Real T_e = T_e_constant;
            if (background_density_is_parsed || background_temperature_is_parsed) {
                amrex::ParticleReal x, y, z;
                GetPosition.AsStored(ip, x, y, z);
                if (m_background_density_is_parsed) n_e = n_e_func(x, y, z, t);
                if (m_background_temperature_is_parsed) T_e = T_e_func(x, y, z, t)*PhysConst::kb;
            }

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

            amrex::Real const vth = sqrt(3._rt*T_e/mass_e);
            amrex::Real const wp = sqrt(n_e*q_e2/(ep0*mass_e));
            amrex::Real const lambdadb = vth/wp;
            amrex::Real const lambdadb3 = lambdadb*lambdadb*lambdadb;
            amrex::Real const loglambda = log((12._rt*pi/Zb)*(n_e*lambdadb3));

            amrex::Real const pi32 = pi*sqrt(pi);
            amrex::Real const q2 = species_charge*species_charge;
            amrex::Real const T32 = T_e*sqrt(T_e);

            amrex::Real alpha = sqrt(2._rt)*n_e*q2*q_e2*sqrt(mass_e)*loglambda/(12._rt*pi32*ep02*species_mass*T32);

            ux[ip] *= exp(-alpha*dt);
            uy[ip] *= exp(-alpha*dt);
            uz[ip] *= exp(-alpha*dt);

        }
        );
}
