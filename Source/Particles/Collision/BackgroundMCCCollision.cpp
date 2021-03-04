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

    // AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_species_names.size() == 2,
    //           "Pair wise Coulomb must have exactly two species.");

    amrex::ParmParse pp(collision_name);

    // default Coulomb log, if < 0, will be computed automatically
    // m_CoulombLog = -1.0_rt;
    // queryWithParser(pp, "CoulombLog", m_CoulombLog);
    pp.query("background_density", m_background_density);

    // query for a list of collision processes
    // collision_types = ['elastic', 'ionization', etc.]
    // build maximum cross-section and collision probability
    // and everything else done in the __init__ function for Warp MCC
    // for testing now I'll just specify a collision probability
    m_total_collision_prob = 0.08356;
    // and use dummy constant cross-section and assume max velocity of 5 keV
    amrex::Real sigma_E = 1e-18;
    m_nu_max = 473299989.9733926;
    m_background_mass = 6.646475849910765e-27;
    m_background_temperature = 300;

}

void
BackgroundMCCCollision::doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc)
{
    WARPX_PROFILE("BackgroundMCCCollision::doCollisions()");

    const amrex::Real dt = WarpX::GetInstance().getdt(0);
    if ( int(std::floor(cur_time/dt)) % m_ndt != 0 ) return;

    auto& species1 = mypc->GetParticleContainerFromName(m_species_names[0]);
    // TODO add functionality to identify secondary species

    m_mass1 = species1.getMass();

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
            doBackgroundCollisionsWithinTile(
                pti, species1
            );
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
    ( WarpXParIter& pti, WarpXParticleContainer& species1)
{
    // get particle count
    const long np = pti.numParticles();

    // get collider properties
    amrex::Real mass1 = m_mass1;

    // get neutral properties
    amrex::Real n_a = m_background_density;
    amrex::Real mass_a = m_background_mass;
    amrex::Real T_a = m_background_temperature;
    amrex::Real vel_std = std::sqrt(PhysConst::kb * T_a / mass_a);

    // get collision parameters
    amrex::Real total_collision_prob = m_total_collision_prob;
    amrex::Real nu_max = m_nu_max;

    // instantiate a random number generator
    amrex::RandomEngine const engine;

    // Get Struct-Of-Array particle data, also called attribs
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();


    // use dummy constant cross-section
    amrex::Real sigma_E = 1e-18;


    amrex::ParallelFor( np,
        [=] AMREX_GPU_DEVICE (long ip)
        {
            // determine if this particle should collide
            if (amrex::Random(engine) < total_collision_prob){

                amrex::Real v_coll, v_coll2, E_coll, nu_i = 0;
                amrex::Real col_select = amrex::Random(engine);
                amrex::ParticleReal ua_x, ua_y, ua_z;
                amrex::ParticleReal uCOM_x, uCOM_y, uCOM_z;

                // FIXME we really should be consistent in using stationary
                // targets for electron collisions.

                // sample random velocity for the target neutral from
                // the neutral temperature
                ua_x = vel_std * amrex::RandomNormal(0, 1.0, engine);
                ua_y = vel_std * amrex::RandomNormal(0, 1.0, engine);
                ua_z = vel_std * amrex::RandomNormal(0, 1.0, engine);

                // calculate the center of momentum velocity
                uCOM_x = (mass1 * ux[ip] + mass_a * ua_x) / (mass1 + mass_a);
                uCOM_y = (mass1 * uy[ip] + mass_a * ua_y) / (mass1 + mass_a);
                uCOM_z = (mass1 * uz[ip] + mass_a * ua_z) / (mass1 + mass_a);

                // calculate collision velocity magnitude and its square
                // for electrons we assume the neutrals are stationary but
                // for ions we have to sample from the neutral temperature
                if (mass_a / mass1 > 1e3){
                    v_coll2 = ux[ip]*ux[ip] + uy[ip]*uy[ip] * uz[ip]*uz[ip];
                }
                else{
                    v_coll2 = (
                        (ux[ip] - ua_x)*(ux[ip] - ua_x)
                        + (uy[ip] - ua_y)*(uy[ip] - ua_y)
                        + (uz[ip] - ua_z)*(uz[ip] - ua_z)
                    );
                }
                v_coll = std::sqrt(v_coll2);

                // calculate particle energy in eV
                // E_coll = 0.5 * mass1 * v_coll2 / PhysConst::q_e;

                // loop through all collision pathways
                // TODO for col_pathway in col_pathways:

                    // get collision cross-section
                    //sigma_E = ...

                    // calculate collision frequency
                    nu_i += n_a * sigma_E * v_coll / nu_max;

                    // check if this collision should be performed and call
                    // the appropriate scattering function
                    if (col_select < nu_i){
                        // if pathway is elastic collision with heavy
                        ElasticScattering(
                            ux[ip], uy[ip], uz[ip], uCOM_x, uCOM_y, uCOM_z,
                            engine
                        );
                    }
            }
        }
    );
}
