/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleCreationFunc.H"

#include "BinaryCollisionUtils.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <string>

ParticleCreationFunc::ParticleCreationFunc (const std::string collision_name,
                                            MultiParticleContainer const * const mypc)
    {
        amrex::ParmParse pp_collision_name(collision_name);

        m_collision_type = BinaryCollisionUtils::get_collision_type(collision_name, mypc);

#ifdef AMREX_USE_GPU
        amrex::Vector<int> host_num_products;
#else
        amrex::Gpu::DeviceVector<int>& host_num_products = m_num_products;
#endif
        if (m_collision_type == CollisionType::ProtonBoronFusion)
            {
                // Proton-Boron fusion only produces alpha particles
                m_num_product_species = 1;
                // Proton-Boron fusion produces 3 alpha particles per fusion reaction
                host_num_products.push_back(3);
            }
#ifdef AMREX_USE_GPU
        m_num_products.resize(m_num_product_species);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, host_num_products.begin(),
                              host_num_products.end(),
                              m_num_products.begin());
        amrex::Gpu::streamSynchronize();
#endif
    }
