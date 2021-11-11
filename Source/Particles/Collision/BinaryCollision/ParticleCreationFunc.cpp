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

        if (m_collision_type == CollisionType::ProtonBoronFusion)
            {
                // Proton-Boron fusion only produces alpha particles
                m_num_product_species = 1;
                // Proton-Boron fusion produces 3 alpha particles per fusion reaction
                m_num_products_host.push_back(3);
#ifndef AMREX_USE_GPU
                // On CPU, the device vector can be filled immediatly
                m_num_products_device.push_back(3);
#endif
            }
#ifdef AMREX_USE_GPU
        m_num_products_device.resize(m_num_product_species);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_num_products_host.begin(),
                              m_num_products_host.end(),
                              m_num_products_device.begin());
        amrex::Gpu::streamSynchronize();
#endif
    }
