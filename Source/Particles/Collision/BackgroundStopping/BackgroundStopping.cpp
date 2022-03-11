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

    std::string background_type_str;
    pp_collision_name.get("background_type", background_type_str);
    if (background_type_str == "electrons") {
        m_background_type = BackgroundStoppingType::ELECTRONS;
    } else if (background_type_str == "ions") {
        m_background_type = BackgroundStoppingType::IONS;
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "background_type must be either electrons or ions");
    }

    amrex::Real background_density;
    std::string background_density_str;
    if (queryWithParser(pp_collision_name, "background_density", background_density)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(background_density > 0,
                 "For background stopping, the background density must be greater than 0");
        m_background_density_parser = makeParser(std::to_string(background_density), {"x", "y", "z", "t"});
    } else if (pp_collision_name.query("background_density(x,y,z,t)", background_density_str)) {
        m_background_density_parser = makeParser(background_density_str, {"x", "y", "z", "t"});
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                 "For background stopping, the background density must be specified.");
    }

    amrex::Real background_temperature;
    std::string background_temperature_str;
    if (queryWithParser(pp_collision_name, "background_temperature", background_temperature)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(background_temperature > 0,
                 "For background stopping, the background temperature must be greater than 0");
        m_background_temperature_parser = makeParser(std::to_string(background_temperature), {"x", "y", "z", "t"});
    } else if (pp_collision_name.query("background_temperature(x,y,z,t)", background_temperature_str)) {
        m_background_temperature_parser = makeParser(background_temperature_str, {"x", "y", "z", "t"});
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                 "For background stopping, the background temperature must be specified.");
    }

    m_background_density_func = m_background_density_parser.compile<4>();
    m_background_temperature_func = m_background_temperature_parser.compile<4>();

    if (m_background_type == BackgroundStoppingType::ELECTRONS) {
        m_background_mass = PhysConst::m_e;
        queryWithParser(pp_collision_name, "background_mass", m_background_mass);
    } else if (m_background_type == BackgroundStoppingType::IONS) {
        getWithParser(pp_collision_name, "background_mass", m_background_mass);
        getWithParser(pp_collision_name, "background_charge_state", m_background_charge_state);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_background_mass > 0,
             "For background stopping, the background mass must be greater than 0");

}

void
BackgroundStopping::doCollisions (amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc)
{
    WARPX_PROFILE("BackgroundStopping::doCollisions()");
    using namespace amrex::literals;

    auto& species = mypc->GetParticleContainerFromName(m_species_names[0]);
    amrex::Real species_mass = species.getMass();
    amrex::Real species_charge = species.getCharge();

    BackgroundStoppingType background_type = m_background_type;

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

            if (background_type == BackgroundStoppingType::ELECTRONS) {
                doBackgroundStoppingOnElectronsWithinTile(pti, dt, cur_time, species_mass, species_charge);
            } else if (background_type == BackgroundStoppingType::IONS) {
                doBackgroundStoppingOnIonsWithinTile(pti, dt, cur_time, species_mass, species_charge);
            }

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
            }
        }

    }
}

void BackgroundStopping::doBackgroundStoppingOnElectronsWithinTile (WarpXParIter& pti, amrex::Real dt, amrex::Real t,
                                                                    amrex::Real species_mass, amrex::Real species_charge)
{
    using namespace amrex::literals;

    // So that CUDA code gets its intrinsic, not the host-only C++ library version
    using std::sqrt, std::abs, std::log, std::exp;

    // get particle count
    const long np = pti.numParticles();

    // get background particle mass
    amrex::Real mass_e = m_background_mass;

    // setup parsers for the background density and temperature
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

            amrex::ParticleReal x, y, z;
            GetPosition.AsStored(ip, x, y, z);
            amrex::Real const n_e = n_e_func(x, y, z, t);
            amrex::Real const T_e = T_e_func(x, y, z, t)*PhysConst::kb;

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

            amrex::Real const alpha = sqrt(2._rt)*n_e*q2*q_e2*sqrt(mass_e)*loglambda/(12._rt*pi32*ep02*species_mass*T32);

            ux[ip] *= exp(-alpha*dt);
            uy[ip] *= exp(-alpha*dt);
            uz[ip] *= exp(-alpha*dt);

        }
        );
}

void BackgroundStopping::doBackgroundStoppingOnIonsWithinTile (WarpXParIter& pti, amrex::Real dt, amrex::Real t,
                                                               amrex::Real species_mass, amrex::Real species_charge)
{
    using namespace amrex::literals;

    // So that CUDA code gets its intrinsic, not the host-only C++ library version
    using std::sqrt, std::abs, std::log, std::exp, std::pow;

    // get particle count
    const long np = pti.numParticles();

    // get background particle mass
    amrex::Real mass_i = m_background_mass;
    amrex::Real charge_state_i = m_background_charge_state;

    // setup parsers for the background density and temperature
    auto n_i_func = m_background_density_func;
    auto T_i_func = m_background_temperature_func;

    // get Struct-Of-Array particle data, also called attribs
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    // May be needed to evaluate the density function
    auto GetPosition = GetParticlePosition(pti);

    amrex::ParallelFor(np,
        [=] AMREX_GPU_HOST_DEVICE (long ip)
        {

            amrex::ParticleReal x, y, z;
            GetPosition.AsStored(ip, x, y, z);
            amrex::Real const n_i = n_i_func(x, y, z, t);
            amrex::Real const T_i = T_i_func(x, y, z, t)*PhysConst::kb;

            // This implements the equation 14.20 from Introduction to Plasma Physics,
            // Goldston and Rutherford, the slowing down of beam ions due to collisions with electrons.
            // The equation is written with energy, W, as dW/dt = -alpha/W**0.5, and integrated to
            // give W(t+dt) = (W(t)**1.5 - 3./2.*alpha*dt)**(2/3)

            amrex::Real constexpr pi = MathConst::pi;
            amrex::Real constexpr q_e = PhysConst::q_e;
            amrex::Real constexpr q_e2 = q_e*q_e;
            amrex::Real constexpr ep0 = PhysConst::ep0;
            amrex::Real constexpr ep02 = ep0*ep0;

            amrex::Real const qi2 = charge_state_i*charge_state_i*q_e2;
            amrex::Real const qb2 = species_charge*species_charge;
            amrex::Real const Zb = abs(species_charge/q_e);

            amrex::Real const vth = sqrt(3._rt*T_i/mass_i);
            amrex::Real const wp = sqrt(n_i*q_e2/(ep0*mass_i));
            amrex::Real const lambdadb = vth/wp;
            amrex::Real const lambdadb3 = lambdadb*lambdadb*lambdadb;
            amrex::Real const loglambda = log((12._rt*pi/Zb)*(n_i*lambdadb3));

            amrex::Real const alpha = sqrt(2._rt)*n_i*qi2*qb2*sqrt(species_mass)*loglambda/(8._rt*pi*ep02*mass_i);

            amrex::Real const W0 = 0.5_rt*species_mass*(ux[ip]*ux[ip] + uy[ip]*uy[ip] + uz[ip]*uz[ip]);
            amrex::Real const f1 = pow(W0, 1.5_rt) - 1.5_rt*alpha*dt;
            // If f1 goes negative, the particle has fully stopped, so set W1 to 0.
            amrex::Real const W1 = pow((f1 > 0._rt ? f1 : 0._rt), 2._rt/3._rt);
            amrex::Real const vscale = (W0 > 0._rt ? std::sqrt(W1/W0) : 0._rt);

            ux[ip] *= vscale;
            uy[ip] *= vscale;
            uz[ip] *= vscale;

        }
        );
}
