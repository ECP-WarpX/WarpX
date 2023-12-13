/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef DEFAULTINITIALIZATION_H_
#define DEFAULTINITIALIZATION_H_

#include <WarpX.H>
#ifdef WARPX_QED
#   include "Particles/ElementaryProcess/QEDInternals/BreitWheelerEngineWrapper.H"
#   include "Particles/ElementaryProcess/QEDInternals/QuantumSyncEngineWrapper.H"
#endif

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <map>
#include <string>

/**
 * \brief This set of initialization policies describes what happens
 * when we need to create a new particle due to an elementary process.
 * For example, when an ionization event creates an electron, these
 * policies control the initial values of the electron's components.
 * These can always be over-written later.
 *
 * The specific meanings are as follows:
 *     Zero         - set the component to zero
 *     One          - set the component to one
 *     RandomExp    - a special flag for the optical depth component used by
 *                    certain QED processes, which gets a random initial value
 *                    extracted from an exponential distribution
 *
 */
enum struct InitializationPolicy {Zero=0, One, RandomExp};

/**
 * \brief This map sets the initialization policy for each particle component
 * used in WarpX.
 */
static std::map<std::string, InitializationPolicy> initialization_policies = {
    {"w",     InitializationPolicy::Zero },
    {"ux",    InitializationPolicy::Zero },
    {"uy",    InitializationPolicy::Zero },
    {"uz",    InitializationPolicy::Zero },
#ifdef WARPX_DIM_RZ
    {"theta", InitializationPolicy::Zero},
#endif

#ifdef WARPX_QED
    {"opticalDepthBW",   InitializationPolicy::RandomExp},
    {"opticalDepthQSR",   InitializationPolicy::RandomExp}
#endif

};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::ParticleReal initializeRealValue (const InitializationPolicy policy, amrex::RandomEngine const& engine) noexcept
{
    switch (policy) {
        case InitializationPolicy::Zero : return 0.0;
        case InitializationPolicy::One  : return 1.0;
        case InitializationPolicy::RandomExp : {
            return -std::log(amrex::Random(engine));
        }
        default : {
            amrex::Abort("Initialization Policy not recognized");
            return 1.0;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int initializeIntValue (const InitializationPolicy policy) noexcept
{
    switch (policy) {
        case InitializationPolicy::Zero : return 0;
        case InitializationPolicy::One  : return 1;
        default : {
            amrex::Abort("Initialization Policy not recognized");
            return 1;
        }
    }
}

namespace ParticleCreation {

    /**
     * \brief Default initialize runtime attributes in a tile. This routine does not initialize the
     * first n_external_attr_real real attributes and the first n_external_attr_int integer
     * attributes, which have been in principle externally set elsewhere.
     *
     * @param[inout] ptile the tile in which attributes are initialized
     * @param[in] n_external_attr_real The number of real attributes that have been externally set.
     * These are NOT initialized by this function.
     * @param[in] n_external_attr_int The number of integer attributes that have been externally set.
     * These are NOT initialized by this function.
     * @param[in] engine the random engine, used in initialization of QED optical depths
     */
template <typename PTile>
void DefaultInitializeRuntimeAttributes (PTile& ptile,
                                         const int n_external_attr_real,
                                         const int n_external_attr_int,
                                         const std::vector<std::string>& user_real_attribs,
                                         const std::vector<std::string>& user_int_attribs,
                                         const std::map<std::string, int>& particle_comps,
                                         const std::map<std::string, int>& particle_icomps,
                                         const std::vector<amrex::Parser*>& user_real_attrib_parser,
                                         const std::vector<amrex::Parser*>& user_int_attrib_parser,
#ifdef WARPX_QED
                                         BreitWheelerEngine* p_bw_engine,
                                         QuantumSynchrotronEngine* p_qs_engine,
#endif
                                         const int ionization_initial_level,
                                         const amrex::RandomEngine& engine)
{
        using namespace amrex::literals;

        const int np = ptile.numParticles();

        // Preparing data needed for user defined attributes
        const auto n_user_real_attribs = static_cast<int>(user_real_attribs.size());
        const auto n_user_int_attribs = static_cast<int>(user_int_attribs.size());
        const auto get_position = GetParticlePosition<PIdx>(ptile);
        const auto soa = ptile.getParticleTileData();
        const amrex::ParticleReal* AMREX_RESTRICT ux = soa.m_rdata[PIdx::ux];
        const amrex::ParticleReal* AMREX_RESTRICT uy = soa.m_rdata[PIdx::uy];
        const amrex::ParticleReal* AMREX_RESTRICT uz = soa.m_rdata[PIdx::uz];
        constexpr int lev = 0;
        const amrex::Real t = WarpX::GetInstance().gett_new(lev);

#ifndef WARPX_QED
        amrex::ignore_unused(engine);
#endif

        // Initialize the last NumRuntimeRealComps() - n_external_attr_real runtime real attributes
        for (int j = PIdx::nattribs + n_external_attr_real; j < ptile.NumRealComps() ; ++j)
        {
            amrex::Vector<amrex::ParticleReal> attr_temp(np, 0.0_prt);
#ifdef WARPX_QED
            // Current runtime comp is quantum synchrotron optical depth
            if (particle_comps.find("opticalDepthQSR") != particle_comps.end() &&
                particle_comps.at("opticalDepthQSR") == j)
            {
                const QuantumSynchrotronGetOpticalDepth quantum_sync_get_opt =
                                                p_qs_engine->build_optical_depth_functor();;
                for (int i = 0; i < np; ++i) {
                    attr_temp[i] = quantum_sync_get_opt(engine);
                }
            }

             // Current runtime comp is Breit-Wheeler optical depth
            if (particle_comps.find("opticalDepthBW") != particle_comps.end() &&
                particle_comps.at("opticalDepthBW") == j)
            {
                const BreitWheelerGetOpticalDepth breit_wheeler_get_opt =
                                                p_bw_engine->build_optical_depth_functor();;
                for (int i = 0; i < np; ++i) {
                    attr_temp[i] = breit_wheeler_get_opt(engine);
                }
            }
#endif

            for (int ia = 0; ia < n_user_real_attribs; ++ia)
            {
                // Current runtime comp is ia-th user defined attribute
                if (particle_comps.find(user_real_attribs[ia]) != particle_comps.end() &&
                    particle_comps.at(user_real_attribs[ia]) == j)
                {
                    amrex::ParticleReal xp, yp, zp;
                    const amrex::ParserExecutor<7> user_real_attrib_parserexec =
                                             user_real_attrib_parser[ia]->compile<7>();
                    for (int i = 0; i < np; ++i) {
                        get_position(i, xp, yp, zp);
                        attr_temp[i] = user_real_attrib_parserexec(xp, yp, zp,
                                                                   ux[i], uy[i], uz[i], t);
                    }
                }
            }

            ptile.push_back_real(j, attr_temp.data(), attr_temp.data() + np);
        }

        // Initialize the last NumRuntimeIntComps() - n_external_attr_int runtime int attributes
        for (int j = n_external_attr_int; j < ptile.NumIntComps() ; ++j)
        {
            amrex::Vector<int> attr_temp(np, 0);

            // Current runtime comp is ionization level
            if (particle_icomps.find("ionizationLevel") != particle_icomps.end() &&
                particle_icomps.at("ionizationLevel") == j)
            {
                for (int i = 0; i < np; ++i) {
                    attr_temp[i] = ionization_initial_level;
                }
            }

            for (int ia = 0; ia < n_user_int_attribs; ++ia)
            {
                // Current runtime comp is ia-th user defined attribute
                if (particle_icomps.find(user_int_attribs[ia]) != particle_icomps.end() &&
                    particle_icomps.at(user_int_attribs[ia]) == j)
                {
                    amrex::ParticleReal xp, yp, zp;
                    const amrex::ParserExecutor<7> user_int_attrib_parserexec =
                                             user_int_attrib_parser[ia]->compile<7>();
                    for (int i = 0; i < np; ++i) {
                        get_position(i, xp, yp, zp);
                        attr_temp[i] = static_cast<int>(
                                user_int_attrib_parserexec(xp, yp, zp, ux[i], uy[i], uz[i], t));
                    }
                }
            }

            ptile.push_back_int(j, attr_temp.data(), attr_temp.data() + np);
        }
}

}

#endif
