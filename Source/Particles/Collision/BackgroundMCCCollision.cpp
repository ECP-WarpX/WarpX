/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCCollision.H"
#include "MCCScattering.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

BackgroundMCCCollision::BackgroundMCCCollision (std::string const collision_name)
    : CollisionBase(collision_name)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() == 1,
                                     "Background MCC must have exactly one species.");

    amrex::ParmParse pp_collision_name(collision_name);

    amrex::Real background_density = 0;
    if (queryWithParser(pp_collision_name, "background_density", background_density)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            (background_density > 0), "The background density must be greater than 0."
        );
        m_background_density_parser = makeParser(std::to_string(background_density), {"x", "y", "z", "t"});
    }
    else {
        std::string background_density_str;
        pp_collision_name.get("background_density(x,y,z,t)", background_density_str);
        m_background_density_parser = makeParser(background_density_str, {"x", "y", "z", "t"});
    }

    amrex::Real background_temperature;
    if (queryWithParser(pp_collision_name, "background_temperature", background_temperature)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            (background_temperature >= 0), "The background temperature must be positive."
        );
        m_background_temperature_parser = makeParser(std::to_string(background_temperature), {"x", "y", "z", "t"});
    }
    else {
        std::string background_temperature_str;
        pp_collision_name.get("background_temperature(x,y,z,t)", background_temperature_str);
        m_background_temperature_parser = makeParser(background_temperature_str, {"x", "y", "z", "t"});
    }

    // compile parsers for background density and temperature
    m_background_density_func = m_background_density_parser.compile<4>();
    m_background_temperature_func = m_background_temperature_parser.compile<4>();

    queryWithParser(pp_collision_name, "max_background_density", m_max_background_density);
    // if the background density is constant we can use that number to calculate
    // the maximum collision probability, if `max_background_density` was not
    // specified
    if (m_max_background_density == 0 && background_density != 0) {
        m_max_background_density = background_density;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (m_max_background_density > 0),
        "The maximum background density must be greater than 0."
    );

    // if the neutral mass is specified use it, but if ionization is
    // included the mass of the secondary species of that interaction
    // will be used. If no neutral mass is specified and ionization is not
    // included the mass of the colliding species will be used
    m_background_mass = -1;
    queryWithParser(pp_collision_name, "background_mass", m_background_mass);

    // query for a list of collision processes
    // these could be elastic, excitation, charge_exchange, back, etc.
    amrex::Vector<std::string> scattering_process_names;
    pp_collision_name.queryarr("scattering_processes", scattering_process_names);

    // create a vector of MCCProcess objects from each scattering
    // process name
    for (const auto& scattering_process : scattering_process_names) {
        std::string kw_cross_section = scattering_process + "_cross_section";
        std::string cross_section_file;
        pp_collision_name.query(kw_cross_section.c_str(), cross_section_file);

        amrex::Real energy = 0.0;
        // if the scattering process is excitation or ionization get the
        // energy associated with that process
        if (scattering_process.find("excitation") != std::string::npos ||
            scattering_process.find("ionization") != std::string::npos) {
            std::string kw_energy = scattering_process + "_energy";
            getWithParser(pp_collision_name, kw_energy.c_str(), energy);
        }

        MCCProcess process(scattering_process, cross_section_file, energy);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(process.type() != MCCProcessType::INVALID,
                                         "Cannot add an unknown MCC process type");

        // if the scattering process is ionization get the secondary species
        // only one ionization process is supported, the vector
        // m_ionization_processes is only used to make it simple to calculate
        // the maximum collision frequency with the same function used for
        // particle conserving processes
        if (process.type() == MCCProcessType::IONIZATION) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!ionization_flag,
                                             "Background MCC only supports a single ionization process");
            ionization_flag = true;

            std::string secondary_species;
            pp_collision_name.get("ionization_species", secondary_species);
            m_species_names.push_back(secondary_species);

            m_ionization_processes.push_back(std::move(process));
        } else {
            m_scattering_processes.push_back(std::move(process));
        }
    }

#ifdef AMREX_USE_GPU
    amrex::Gpu::HostVector<MCCProcess::Executor> h_scattering_processes_exe;
    amrex::Gpu::HostVector<MCCProcess::Executor> h_ionization_processes_exe;
    for (auto const& p : m_scattering_processes) {
        h_scattering_processes_exe.push_back(p.executor());
    }
    for (auto const& p : m_ionization_processes) {
        h_ionization_processes_exe.push_back(p.executor());
    }
    m_scattering_processes_exe.resize(h_scattering_processes_exe.size());
    m_ionization_processes_exe.resize(h_ionization_processes_exe.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_scattering_processes_exe.begin(),
                          h_scattering_processes_exe.end(), m_scattering_processes_exe.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_ionization_processes_exe.begin(),
                          h_ionization_processes_exe.end(), m_ionization_processes_exe.begin());
    amrex::Gpu::streamSynchronize();
#else
    for (auto const& p : m_scattering_processes) {
        m_scattering_processes_exe.push_back(p.executor());
    }
    for (auto const& p : m_ionization_processes) {
        m_ionization_processes_exe.push_back(p.executor());
    }
#endif
}

/** Calculate the maximum collision frequency using a fixed energy grid that
 *  ranges from 1e-4 to 5000 eV in 0.2 eV increments
 */
amrex::Real
BackgroundMCCCollision::get_nu_max(amrex::Vector<MCCProcess> const& mcc_processes)
{
    using namespace amrex::literals;
    amrex::Real nu, nu_max = 0.0;

    for (amrex::Real E = 1e-4_rt; E < 5000._rt; E+=0.2_rt) {
        amrex::Real sigma_E = 0.0;

        // loop through all collision pathways
        for (const auto &scattering_process : mcc_processes) {
            // get collision cross-section
            sigma_E += scattering_process.getCrossSection(E);
        }

        // calculate collision frequency
        nu = (
              m_max_background_density
              * std::sqrt(2.0_rt / m_mass1 * PhysConst::q_e)
              * sigma_E * std::sqrt(E)
              );
        if (nu > nu_max) {
            nu_max = nu;
        }
    }
    return nu_max;
}

void
BackgroundMCCCollision::doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc)
{
    WARPX_PROFILE("BackgroundMCCCollision::doCollisions()");
    using namespace amrex::literals;

    const amrex::Real dt = WarpX::GetInstance().getdt(0);
    if ( int(std::floor(cur_time/dt)) % m_ndt != 0 ) return;

    auto& species1 = mypc->GetParticleContainerFromName(m_species_names[0]);
    // this is a very ugly hack to have species2 be a reference and be
    // defined in the scope of doCollisions
    auto& species2 = (
                      (m_species_names.size() == 2) ?
                      mypc->GetParticleContainerFromName(m_species_names[1]) :
                      mypc->GetParticleContainerFromName(m_species_names[0])
                      );

    if (!init_flag) {
        m_mass1 = species1.getMass();

        // calculate maximum collision frequency without ionization
        m_nu_max = get_nu_max(m_scattering_processes);

        // calculate total collision probability
        auto coll_n = m_nu_max * m_ndt * dt;
        m_total_collision_prob = 1.0_rt - std::exp(-coll_n);

        // dt has to be small enough that a linear expansion of the collision
        // probability is sufficiently accurately, otherwise the MCC results
        // will be very heavily affected by small changes in the timestep
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(coll_n < 0.1_rt,
            "dt is too large to ensure accurate MCC results"
        );

        if (ionization_flag) {
            // calculate maximum collision frequency for ionization
            m_nu_max_ioniz = get_nu_max(m_ionization_processes);

            // calculate total ionization probability
            auto coll_n_ioniz = m_nu_max_ioniz * m_ndt * dt;
            m_total_collision_prob_ioniz = 1.0_rt - std::exp(-coll_n_ioniz);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(coll_n_ioniz < 0.1_rt,
                "dt is too large to ensure accurate MCC results"
            );

            // if an ionization process is included the secondary species mass
            // is taken as the background mass
            m_background_mass = species2.getMass();
        }
        // if no neutral species mass was specified and ionization is not
        // included assume that the collisions will be with neutrals of the
        // same mass as the colliding species (like with ion-neutral collisions)
        else if (m_background_mass == -1) {
            m_background_mass = species1.getMass();
        }

        amrex::Print() <<
            "Setting up collisions for " << m_species_names[0] << " with total "
            "collision probability: " <<
            m_total_collision_prob + m_total_collision_prob_ioniz << "\n";

        init_flag = true;
    }

    // Loop over refinement levels
    auto const flvl = species1.finestLevel();
    for (int lev = 0; lev <= flvl; ++lev) {

        auto cost = WarpX::getCosts(lev);

        // firstly loop over particles box by box and do all particle conserving
        // scattering
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(species1, lev); pti.isValid(); ++pti) {
            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            amrex::Real wt = amrex::second();

            doBackgroundCollisionsWithinTile(pti, cur_time);

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add( &(*cost)[pti.index()], wt);
            }
        }

        // secondly perform ionization through the SmartCopyFactory if needed
        if (ionization_flag) {
            doBackgroundIonization(lev, cost, species1, species2, cur_time);
        }
    }
}


void BackgroundMCCCollision::doBackgroundCollisionsWithinTile
( WarpXParIter& pti, amrex::Real t )
{
    using namespace amrex::literals;

    // So that CUDA code gets its intrinsic, not the host-only C++ library version
    using std::sqrt;

    // get particle count
    const long np = pti.numParticles();

    // get collider properties
    amrex::Real mass1 = m_mass1;

    // get neutral mass
    amrex::Real mass_a = m_background_mass;
    // get parsers for the background density and temperature
    auto n_a_func = m_background_density_func;
    auto T_a_func = m_background_temperature_func;

    // get collision parameters
    auto scattering_processes = m_scattering_processes_exe.data();
    int const process_count   = m_scattering_processes_exe.size();

    amrex::Real total_collision_prob = m_total_collision_prob;
    amrex::Real nu_max = m_nu_max;

    // we need particle positions in order to calculate the local density
    // and temperature
    auto GetPosition = GetParticlePosition(pti);

    // get Struct-Of-Array particle data, also called attribs
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    amrex::ParallelForRNG(np,
                          [=] AMREX_GPU_HOST_DEVICE (long ip, amrex::RandomEngine const& engine)
                          {
                              // determine if this particle should collide
                              if (amrex::Random(engine) > total_collision_prob) return;

                              amrex::ParticleReal x, y, z;
                              GetPosition.AsStored(ip, x, y, z);

                              amrex::Real n_a = n_a_func(x, y, z, t);
                              amrex::Real T_a = T_a_func(x, y, z, t);

                              amrex::Real vel_std, v_coll, v_coll2, E_coll, sigma_E, nu_i = 0;
                              amrex::Real col_select = amrex::Random(engine);
                              amrex::ParticleReal ua_x, ua_y, ua_z;
                              amrex::ParticleReal uCOM_x, uCOM_y, uCOM_z;

                              // get velocities of gas particles from a Maxwellian distribution
                              vel_std = sqrt(PhysConst::kb * T_a / mass_a);
                              ua_x = vel_std * amrex::RandomNormal(0_rt, 1.0_rt, engine);
                              ua_y = vel_std * amrex::RandomNormal(0_rt, 1.0_rt, engine);
                              ua_z = vel_std * amrex::RandomNormal(0_rt, 1.0_rt, engine);

                              // calculate the center of momentum velocity
                              uCOM_x = (mass1 * ux[ip] + mass_a * ua_x) / (mass1 + mass_a);
                              uCOM_y = (mass1 * uy[ip] + mass_a * ua_y) / (mass1 + mass_a);
                              uCOM_z = (mass1 * uz[ip] + mass_a * ua_z) / (mass1 + mass_a);

                              // calculate relative velocity of collision and collision energy if
                              // the colliding particle is an ion. For electron collisions we
                              // cannot use the relative velocity since that allows the
                              // possibility where the electron kinetic energy in the lab frame
                              // is insufficient to cause excitation but not in the COM frame -
                              // for energy to balance this situation requires the neutral to
                              // lose energy during the collision which we don't currently
                              // account for.
                              if (mass_a / mass1 > 1e3) {
                                  v_coll2 = ux[ip]*ux[ip] + uy[ip]*uy[ip] + uz[ip]*uz[ip];
                                  E_coll = 0.5_rt * mass1 * v_coll2 / PhysConst::q_e;
                              }
                              else {
                                  v_coll2 = (
                                             (ux[ip] - ua_x)*(ux[ip] - ua_x)
                                             + (uy[ip] - ua_y)*(uy[ip] - ua_y)
                                             + (uz[ip] - ua_z)*(uz[ip] - ua_z)
                                             );
                                  E_coll = (
                                            0.5_rt * mass1 * mass_a / (mass1 + mass_a) * v_coll2
                                            / PhysConst::q_e
                                            );
                              }
                              v_coll = sqrt(v_coll2);

                              // loop through all collision pathways
                              for (int i = 0; i < process_count; i++) {
                                  auto const& scattering_process = *(scattering_processes + i);

                                  // get collision cross-section
                                  sigma_E = scattering_process.getCrossSection(E_coll);

                                  // calculate normalized collision frequency
                                  nu_i += n_a * sigma_E * v_coll / nu_max;

                                  // check if this collision should be performed and call
                                  // the appropriate scattering function
                                  if (col_select > nu_i) continue;

                                  if (scattering_process.m_type == MCCProcessType::ELASTIC) {
                                      ElasticScattering(
                                                        ux[ip], uy[ip], uz[ip], uCOM_x, uCOM_y, uCOM_z, engine
                                                        );
                                  }
                                  else if (scattering_process.m_type == MCCProcessType::BACK) {
                                      BackScattering(
                                                     ux[ip], uy[ip], uz[ip], uCOM_x, uCOM_y, uCOM_z
                                                     );
                                  }
                                  else if (scattering_process.m_type == MCCProcessType::CHARGE_EXCHANGE) {
                                      ChargeExchange(ux[ip], uy[ip], uz[ip], ua_x, ua_y, ua_z);
                                  }
                                  else if (scattering_process.m_type == MCCProcessType::EXCITATION) {
                                      // get the new velocity magnitude
                                      amrex::Real vp = sqrt(
                                                            2.0_rt / mass1 * PhysConst::q_e
                                                            * (E_coll - scattering_process.m_energy_penalty)
                                                            );
                                      ParticleUtils::RandomizeVelocity(ux[ip], uy[ip], uz[ip], vp, engine);
                                  }
                                  break;
                              }
                          }
                          );
}


void BackgroundMCCCollision::doBackgroundIonization
( int lev, amrex::LayoutData<amrex::Real>* cost,
  WarpXParticleContainer& species1, WarpXParticleContainer& species2, amrex::Real t)
{
    WARPX_PROFILE("BackgroundMCCCollision::doBackgroundIonization()");

    SmartCopyFactory copy_factory_elec(species1, species1);
    SmartCopyFactory copy_factory_ion(species1, species2);
    const auto CopyElec = copy_factory_elec.getSmartCopy();
    const auto CopyIon = copy_factory_ion.getSmartCopy();

    const auto Filter = ImpactIonizationFilterFunc(
                                                   m_ionization_processes[0],
                                                   m_mass1, m_total_collision_prob_ioniz,
                                                   m_nu_max_ioniz, m_background_density_func, t
                                                   );

    amrex::Real sqrt_kb_m = std::sqrt(PhysConst::kb / m_background_mass);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (WarpXParIter pti(species1, lev); pti.isValid(); ++pti) {

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        amrex::Real wt = amrex::second();

        auto& elec_tile = species1.ParticlesAt(lev, pti);
        auto& ion_tile = species2.ParticlesAt(lev, pti);

        const auto np_elec = elec_tile.numParticles();
        const auto np_ion = ion_tile.numParticles();

        auto Transform = ImpactIonizationTransformFunc(
                                                       m_ionization_processes[0].getEnergyPenalty(),
                                                       m_mass1, sqrt_kb_m, m_background_temperature_func, t
                                                       );

        const auto num_added = filterCopyTransformParticles<1>(
                                                               elec_tile, ion_tile, elec_tile, np_elec, np_ion,
                                                               Filter, CopyElec, CopyIon, Transform
                                                               );

        setNewParticleIDs(elec_tile, np_elec, num_added);
        setNewParticleIDs(ion_tile, np_ion, num_added);

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[pti.index()], wt);
        }
    }
}
