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

#include <EmbeddedBoundary/LevelSet.H>

template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& level_set, int lev, F&& f)
{
    scrapeParticles(pc, level_set, lev, lev, std::forward<F>(f));
}

template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& level_set, F&& f)
{
    scrapeParticles(pc, level_set, 0, pc.finestLevel(), std::forward<F>(f));
}

template <class PC, class F, std::enable_if_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
void
scrapeParticles (PC& pc, const amrex::Vector<const amrex::MultiFab*>& level_set,
                 int lev_min, int lev_max, F&& f)
{
    BL_PROFILE("scrapeParticles");

    using ParIter = typename PC::ParIterType;
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        const auto& ls = *level_set[lev];
        const auto& geom = pc.Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& tile = pti.GetParticleTile();
            auto ptd = tile.getParticleTileData();
            const auto np = tile.numParticles();
            auto pstruct = tile.GetArrayOfStructs()().data();
            auto phi = ls[pti].array();  // signed distance function
            amrex::ParallelForRNG( np,
            [=] AMREX_GPU_DEVICE (const int ip, amrex::RandomEngine const& engine) noexcept
            {
                auto& p = pstruct[ip];

                int i, j, k;
                amrex::Real W[3][2];

                LevelSet::compute_weights(p, i, j, k, W, plo, dxi);
                amrex::Real phi_value  = LevelSet::interp_ls(p, i, j, k, W, phi, plo);
                amrex::RealVect normal = LevelSet::interp_normal(p, i, j, k, W, phi, plo, dxi);

                // the closest point on the surface to pos is pos - grad phi(pos) * phi(pos)
                amrex::RealVect pos;
                pos[0] = p.pos(0) - normal[0]*phi_value;
                pos[1] = p.pos(1) - normal[1]*phi_value;
                pos[2] = p.pos(2) - normal[2]*phi_value;

                LevelSet::normalize(normal);

                if (phi_value < 0.0) {
                    f(ptd, ip, pos, normal, engine);
                }
            });
        }
    }
}

#endif
