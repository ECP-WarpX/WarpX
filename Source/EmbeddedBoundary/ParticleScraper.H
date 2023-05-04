/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLESCRAPER_H_
#define PARTICLESCRAPER_H_

#include "EmbeddedBoundary/DistanceToEB.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <ablastr/particles/NodalFieldGather.H>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particle.H>
#include <AMReX_RandomEngine.H>
#include <AMReX_REAL.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Vector.H>

#include <type_traits>
#include <utility>



/**
 * \brief Interact particles with the embedded boundary walls.
 *
 *  This function detects which particles have entered into the region
 *  covered by the embedded boundaries and applies an operation on those
 *  that have. Boundary collision detection is performed using a signed
 *  distance function, which is generated automatically when WarpX is
 *  compiled with EB support.
 *
 *  The operation to be performed is specified by the callable function
 *  passed in to this function as an argument. This function can access the
 *  position at which the particle hit the boundary, and also the associated
 *  normal vector. Particles can be `absorbed` by setting their ids to negative
 *  to flag them for removal. Likewise, the can be reflected back into the domain
 *  by modifying their data appropriately and leaving their ids alone.
 *
 *  This version operates only at the specified level.
 *
 * \tparam pc a type of amrex ParticleContainer
 * \tparam F a callable type, e.g. a lambda function or functor
 *
 * \param pc the particle container to test for boundary interactions.
 * \param distance_to_eb a set of MultiFabs that store the signed distance function
 * \param lev the mesh refinement level to work on.
 * \param f the callable that defines what to do when a particle hits the boundary.
 *
 *        The form of the callable should model:
 *        template <typename PData>
 *        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 *        void operator() (const PData& ptd, int i,
 *                         const amrex::RealVect& pos, const amrex::RealVect& normal,
 *                         amrex::RandomEngine const& engine);
 *
 *        where ptd is the particle tile, i the index of the particle operated on,
 *        pos and normal the location of the collision and the boundary normal vector.
 *        engine is for random number generation, if needed.
 */
template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& distance_to_eb, int lev, F&& f)
{
    scrapeParticles(pc, distance_to_eb, lev, lev, std::forward<F>(f));
}

/**
 * \brief Interact particles with the embedded boundary walls.
 *
 *  This function detects which particles have entered into the region
 *  covered by the embedded boundaries and applies an operation on those
 *  that have. Boundary collision detection is performed using a signed
 *  distance function, which is generated automatically when WarpX is
 *  compiled with EB support.
 *
 *  The operation to be performed is specified by the callable function
 *  passed in to this function as an argument. This function can access the
 *  position at which the particle hit the boundary, and also the associated
 *  normal vector. Particles can be `absorbed` by setting their ids to negative
 *  to flag them for removal. Likewise, the can be reflected back into the domain
 *  by modifying their data appropriately and leaving their ids alone.
 *
 *  This version operates over all the levels in the pc.
 *
 * \tparam pc a type of amrex ParticleContainer
 * \tparam F a callable type, e.g. a lambda function or functor
 *
 * \param pc the particle container to test for boundary interactions.
 * \param distance_to_eb a set of MultiFabs that store the signed distance function
 * \param f the callable that defines what to do when a particle hits the boundary.
 *
 *        The form of the callable should model:
 *        template <typename PData>
 *        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 *        void operator() (const PData& ptd, int i,
 *                         const amrex::RealVect& pos, const amrex::RealVect& normal,
 *                         amrex::RandomEngine const& engine);
 *
 *        where ptd is the particle tile, i the index of the particle operated on,
 *        pos and normal the location of the collision and the boundary normal vector.
 *        engine is for random number generation, if needed.
 */
template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& distance_to_eb, F&& f)
{
    scrapeParticles(pc, distance_to_eb, 0, pc.finestLevel(), std::forward<F>(f));
}

/**
 * \brief Interact particles with the embedded boundary walls.
 *
 *  This function detects which particles have entered into the region
 *  covered by the embedded boundaries and applies an operation on those
 *  that have. Boundary collision detection is performed using a signed
 *  distance function, which is generated automatically when WarpX is
 *  compiled with EB support.
 *
 *  The operation to be performed is specified by the callable function
 *  passed in to this function as an argument. This function can access the
 *  position at which the particle hit the boundary, and also the associated
 *  normal vector. Particles can be `absorbed` by setting their ids to negative
 *  to flag them for removal. Likewise, the can be reflected back into the domain
 *  by modifying their data appropriately and leaving their ids alone.
 *
 *  This version operates only at the specified levels.
 *
 * \tparam pc a type of amrex ParticleContainer
 * \tparam F a callable type, e.g. a lambda function or functor
 *
 * \param pc the particle container to test for boundary interactions.
 * \param distance_to_eb a set of MultiFabs that store the signed distance function
 * \param lev_min the minimum mesh refinement level to work on.
 * \param lev_max the maximum mesh refinement level to work on.
 * \param f the callable that defines what to do when a particle hits the boundary.
 *
 *        The form of the callable should model:
 *        template <typename PData>
 *        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 *        void operator() (const PData& ptd, int i,
 *                         const amrex::RealVect& pos, const amrex::RealVect& normal,
 *                         amrex::RandomEngine const& engine);
 *
 *        where ptd is the particle tile, i the index of the particle operated on,
 *        pos and normal the location of the collision and the boundary normal vector.
 *        engine is for random number generation, if needed.
 */
template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& distance_to_eb,
                 int lev_min, int lev_max, F&& f)
{
    BL_PROFILE("scrapeParticles");

    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        const auto plo = pc.Geom(lev).ProbLoArray();
        const auto dxi = pc.Geom(lev).InvCellSizeArray();
        for(WarpXParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            const auto getPosition = GetParticlePosition(pti);
            auto& tile = pti.GetParticleTile();
            auto ptd = tile.getParticleTileData();
            const auto np = tile.numParticles();
            amrex::Particle<0,0> * const particles = tile.GetArrayOfStructs()().data();
            auto phi = (*distance_to_eb[lev])[pti].array();  // signed distance function
            amrex::ParallelForRNG( np,
            [=] AMREX_GPU_DEVICE (const int ip, amrex::RandomEngine const& engine) noexcept
            {
                // skip particles that are already flagged for removal
                if (particles[ip].id() < 0) return;

                amrex::ParticleReal xp, yp, zp;
                getPosition(ip, xp, yp, zp);

                int i, j, k;
                amrex::Real W[AMREX_SPACEDIM][2];
                ablastr::particles::compute_weights_nodal(xp, yp, zp, plo, dxi, i, j, k, W);

                amrex::Real phi_value  = ablastr::particles::interp_field_nodal(i, j, k, W, phi);

                if (phi_value < 0.0)
                {
                    amrex::RealVect normal = DistanceToEB::interp_normal(i, j, k, W, phi, dxi);

                    // the closest point on the surface to pos is pos - grad phi(pos) * phi(pos)
                    amrex::RealVect pos;
#if (defined WARPX_DIM_3D)
                    pos[0] = xp - normal[0]*phi_value;
                    pos[1] = yp - normal[1]*phi_value;
                    pos[2] = zp - normal[2]*phi_value;
#elif (defined WARPX_DIM_XZ)
                    pos[0] = xp - normal[0]*phi_value;
                    pos[1] = zp - normal[1]*phi_value;
#elif (defined WARPX_DIM_RZ)
                    pos[0] = std::sqrt(xp*xp + yp*yp) - normal[0]*phi_value;
                    pos[1] = zp - normal[1]*phi_value;
#elif (defined WARPX_DIM_1D_Z)
                    pos[0] = zp - normal[0]*phi_value;
#endif
                    DistanceToEB::normalize(normal);

                    f(ptd, ip, pos, normal, engine);
                }
            });
        }
    }
}

#endif
