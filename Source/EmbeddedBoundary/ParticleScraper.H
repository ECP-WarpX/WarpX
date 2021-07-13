/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLESCRAPER_H_
#define PARTICLESCRAPER_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#include <EmbeddedBoundary/DistanceToEB.H>

template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& distance_to_eb, int lev, F&& f)
{
    scrapeParticles(pc, distance_to_eb, lev, lev, std::forward<F>(f));
}

template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& distance_to_eb, F&& f)
{
    scrapeParticles(pc, distance_to_eb, 0, pc.finestLevel(), std::forward<F>(f));
}

template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& distance_to_eb,
                 int lev_min, int lev_max, F&& f)
{
    BL_PROFILE("scrapeParticles");

    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        const auto dxi = pc.Geom(lev).InvCellSizeArray();
        for(WarpXParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            const auto getPosition = GetParticlePosition(pti);
            auto& tile = pti.GetParticleTile();
            auto ptd = tile.getParticleTileData();
            const auto np = tile.numParticles();
            auto phi = (*distance_to_eb[lev])[pti].array();  // signed distance function
            amrex::ParallelForRNG( np,
            [=] AMREX_GPU_DEVICE (const int ip, amrex::RandomEngine const& engine) noexcept
            {
                amrex::ParticleReal xp, yp, zp;
                getPosition(ip, xp, yp, zp);

                int i, j, k;
                amrex::Real W[AMREX_SPACEDIM][2];
                DistanceToEB::compute_weights(xp, yp, zp, i, j, k, W);

                amrex::Real phi_value  = DistanceToEB::interp_distance(i, j, k, W, phi);
                amrex::RealVect normal = DistanceToEB::interp_normal(i, j, k, W, phi, dxi);

                // the closest point on the surface to pos is pos - grad phi(pos) * phi(pos)
                amrex::RealVect pos;
                pos[0] = xp - normal[0]*phi_value;
                pos[1] = yp - normal[1]*phi_value;
                pos[2] = zp - normal[2]*phi_value;

                DistanceToEB::normalize(normal);

                if (phi_value < 0.0) {
                    f(ptd, ip, pos, normal, engine);
                }
            });
        }
    }
}

#endif
