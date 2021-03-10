/* Copyright 2021 Roelof Groenewald
 *
 * This file is part of WarpX.
 *
 * License: ????
 */
#include "BackgroundMCCCollision.H"
#include "MCCScattering.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

// using namespace amrex::literals;

BackgroundMCCCollision::BackgroundMCCCollision (std::string const collision_name)
    : CollisionBase(collision_name)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() == 1,
               "Background MCC must have exactly one species.");

    amrex::ParmParse pp(collision_name);

    pp.query("background_density", m_background_density);
    pp.query("background_temperature", m_background_temperature);

    // if the neutral mass is specified use it, but if ionization is
    // included the mass of the secondary species of that interaction
    // will be used. If no neutral mass is specified and ionization is not
    // included the mass of the colliding species will be used
    m_background_mass = -1;
    pp.query("background_mass", m_background_mass);

    // query for a list of collision processes
    // these could be ionization, elastic, charge_exchange, back, etc.
    amrex::Vector<std::string> scattering_process_names;
    pp.queryarr("scattering_processes", scattering_process_names);

    // create a vector of MCCProcess objects from each scattering
    // process name
    for (auto scattering_process : scattering_process_names) {
        std::string temp = scattering_process;
        temp.append("_cross_section");
        char kw[temp.length() + 1];
        std::strcpy(kw, temp.c_str());
        std::string cross_section_file;
        pp.query(kw, cross_section_file);

        amrex::Real energy = 0.0;
        // if the scattering process is excitation or ionization get the
        // energy associated with that process
        if (scattering_process.find("excitation") != std::string::npos ||
            scattering_process.find("ionization") != std::string::npos)
        {
            std::string temp = scattering_process;
            temp.append("_energy");
            char kw[temp.length() + 1];
            std::strcpy(kw, temp.c_str());
            pp.get(kw, energy);
        }

        // if the scattering process is ionization get the secondary species
        if (scattering_process == "ionization")
        {
            std::string secondary_species;
            pp.get("ionization_species", secondary_species);
            m_species_names.push_back(secondary_species);
        }

        m_scattering_processes.push_back(
            MCCProcess(scattering_process, cross_section_file, energy)
        );
    }
}

/** Calculate the maximum collision frequency using a fixed energy grid that
 *  ranges from 1e-4 to 5000 eV in 0.2 eV increments
 */
amrex::Real
BackgroundMCCCollision::get_nu_max()
{
    amrex::Real nu, nu_max = 0.0;

    for (double E = 1e-4; E < 5000; E+=0.2)
    {
        amrex::Real sigma_E = 0.0;

        // loop through all collision pathways
        for (auto scattering_process : m_scattering_processes){
            // get collision cross-section
            sigma_E += scattering_process.getCrossSection(E);
        }

        // calculate collision frequency
        nu = (
            m_background_density * std::sqrt(2.0 / m_mass1 * PhysConst::q_e)
            * sigma_E * std::sqrt(E)
        );
        if (nu > nu_max)
        {
            nu_max = nu;
        }
    }
    return nu_max;
}

void
BackgroundMCCCollision::doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc)
{
    WARPX_PROFILE("BackgroundMCCCollision::doCollisions()");

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

    if (!init_flag)
    {
        m_mass1 = species1.getMass();

        // calculate maximum collision frequency
        m_nu_max = get_nu_max();

        // calculate total collision probability
        auto coll_n = m_nu_max * dt;
        m_total_collision_prob = 1.0 - std::exp(-coll_n);

        // dt has to be small enough that a linear expansion of the collision
        // probability is sufficiently accurately, otherwise the MCC results
        // will be very heavily affected by small changes in the timestep
        if (coll_n > 0.1){
            amrex::Print() <<
                "dt is too large to ensure accurate MCC results for "
                << m_species_names[0] << " collisions since nu_max*dt = "
                << coll_n << " > 0.1\n";
            amrex::Abort();
        }
        amrex::Print() <<
            "Setting up collisions for " << m_species_names[0] << " with total "
            "collision probability: " << m_total_collision_prob << "\n";

        // if an ionization process is included the secondary species mass is
        // taken as the background mass
        if (m_species_names.size() == 2){
            m_background_mass = species2.getMass();
        }
        // if no neutral species mass was specified and ionization is not
        // included assume that the collisions will be with neutrals of the
        // same mass as the colliding species (like with ion-neutral collisions)
        else if (m_background_mass == -1){
            m_background_mass = species1.getMass();
        }

        init_flag = true;
    }

    // Enable tiling
    // amrex::MFItInfo info;
    // if (amrex::Gpu::notInLaunchRegion()) info.EnableTiling(species1.tile_size);

    // Loop over refinement levels
    for (int lev = 0; lev <= species1.finestLevel(); ++lev){

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Loop over particles, box by box
        for (WarpXParIter pti(species1, lev); pti.isValid(); ++pti)
        {
            doBackgroundCollisionsWithinTile(pti, species1, species2);
        }

        /*
        // Loop over all grids/tiles at this level
#ifdef AMREX_USE_OMP
        info.SetDynamic(true);
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi = species1.MakeMFIter(lev, info); mfi.isValid(); ++mfi){
            doBackgroundCollisionsWithinTile( lev, mfi, species1 );
        }
        */
    }
}


// using namespace ParticleUtils;

/** Perform all MCC collisions within a tile
 *
 * @param pti particle iterator
 * @param species1/2 pointer to species container used to inject new particles
          from ionization events
 *
 */
void BackgroundMCCCollision::doBackgroundCollisionsWithinTile
    ( WarpXParIter& pti, WarpXParticleContainer& species1,
      WarpXParticleContainer& species2)
{
    // get particle count
    const long np = pti.numParticles();

    // get collider properties
    amrex::Real mass1 = m_mass1;

    // get neutral properties
    amrex::Real n_a = m_background_density;
    amrex::Real T_a = m_background_temperature;
    amrex::Real mass_a = m_background_mass;
    amrex::Real vel_std = std::sqrt(PhysConst::kb * T_a / mass_a);

    // get collision parameters
    auto scattering_processes = m_scattering_processes;
    amrex::Real total_collision_prob = m_total_collision_prob;
    amrex::Real nu_max = m_nu_max;

    // Get Struct-Of-Array particle data, also called attribs
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    amrex::ParallelForRNG( np,
        [=] AMREX_GPU_DEVICE (long ip, amrex::RandomEngine const& engine)
        {
            // determine if this particle should collide
            if (amrex::Random(engine) > total_collision_prob) return;

            amrex::Real v_coll, v_coll2, E_coll, sigma_E, nu_i = 0;
            amrex::Real col_select = amrex::Random(engine);
            amrex::ParticleReal ua_x, ua_y, ua_z;
            amrex::ParticleReal uCOM_x, uCOM_y, uCOM_z;

            // FIXME we really should be consistent in using stationary
            // targets for electron collisions.

            // get velocities of gas particles from a Maxwellian distribution
            ua_x = vel_std * amrex::RandomNormal(0, 1.0, engine);
            ua_y = vel_std * amrex::RandomNormal(0, 1.0, engine);
            ua_z = vel_std * amrex::RandomNormal(0, 1.0, engine);

            // calculate the center of momentum velocity
            uCOM_x = (mass1 * ux[ip] + mass_a * ua_x) / (mass1 + mass_a);
            uCOM_y = (mass1 * uy[ip] + mass_a * ua_y) / (mass1 + mass_a);
            uCOM_z = (mass1 * uz[ip] + mass_a * ua_z) / (mass1 + mass_a);

            // calculate relative velocity of collision and collision energy if
            // the colliding particle is an ion. For electron collisions we
            // cannot use the relative velocity since that allows the
            // possibility where the electron kinetic energy in the lab frame
            // is insufficient to cause ionization but not in the COM frame -
            // for energy to balance this situation requires the neutral to
            // lose energy during the collision which we don't currently
            // account for.
            if (mass_a / mass1 > 1e3){
                v_coll2 = ux[ip]*ux[ip] + uy[ip]*uy[ip] * uz[ip]*uz[ip];
                E_coll = 0.5 * mass1 * v_coll2 / PhysConst::q_e;
            }
            else{
                v_coll2 = (
                    (ux[ip] - ua_x)*(ux[ip] - ua_x)
                    + (uy[ip] - ua_y)*(uy[ip] - ua_y)
                    + (uz[ip] - ua_z)*(uz[ip] - ua_z)
                );
                E_coll = (
                    0.5 * mass1 * mass_a / (mass1 + mass_a) * v_coll2
                    / PhysConst::q_e
                );
            }
            v_coll = std::sqrt(v_coll2);

            // loop through all collision pathways
            for (auto scattering_process : scattering_processes){

                // get collision cross-section
                sigma_E = scattering_process.getCrossSection(E_coll);

                // calculate collision frequency
                nu_i += n_a * sigma_E * v_coll / nu_max;

                // check if this collision should be performed and call
                // the appropriate scattering function
                if (col_select > nu_i) continue;
                if (scattering_process.name == "elastic"){
                    ElasticScattering(
                        ux[ip], uy[ip], uz[ip], uCOM_x, uCOM_y, uCOM_z, engine
                    );
                }
                else if (scattering_process.name == "back"){
                    BackScattering(
                        ux[ip], uy[ip], uz[ip], uCOM_x, uCOM_y, uCOM_z
                    );
                }
                else if (scattering_process.name == "charge_exchange"){
                    ChargeExchange(ux[ip], uy[ip], uz[ip], ua_x, ua_y, ua_z);
                }
                else if (scattering_process.name.find("excitation")
                        != std::string::npos)
                {
                    // get the new velocity magnitude
                    amrex::Real vp = std::sqrt(
                        2.0 / mass1 * PhysConst::q_e
                        * (E_coll - scattering_process.energy_penalty)
                    );
                    Excitation(ux[ip], uy[ip], uz[ip], vp, engine);
                }
                else if (scattering_process.name == "ionization"){
                    // Ionization(
                    //    ux[ip], uy[ip], uz[ip], ua_x, ua_y, ua_z
                    //);
                    amrex::Print() << "Not implemented yet.";
                }
                break;
            }
        }
    );
}
