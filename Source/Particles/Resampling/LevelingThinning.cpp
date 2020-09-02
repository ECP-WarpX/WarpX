/* Copyright 2019-2020 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "LevelingThinning.H"
#include "Particles/Collision/CollisionType.H"

#include <AMReX_Particles.H>

LevelingThinning::LevelingThinning ()
{
    using namespace amrex::literals;

    amrex::ParmParse pp("resampling_algorithm");
    pp.query("target_ratio", m_target_ratio);
    if (m_target_ratio <= 1._rt)
    {
        amrex::Warning("WARNING: target ratio for leveling thinning is smaller or equal to one."
                       " It is possible that resampling will not remove any particle");
    }
}

LevelingThinning::~LevelingThinning () {}

void LevelingThinning::operator() (WarpXParIter& pti, const int lev, WarpXParticleContainer *pc) const
{
    using namespace amrex::literals;

    auto& ptile = pc->ParticlesAt(lev, pti);

    auto bins = findParticles::findParticlesInEachCell(lev, pti, ptile);

    const int n_cells = bins.numBins();

    auto& soa = ptile.GetStructOfArrays();

    amrex::ParticleReal * const AMREX_RESTRICT w = soa.GetRealData(PIdx::w).data();

    const auto particle_ptr = ptile.GetArrayOfStructs()().data();

    auto indices = bins.permutationPtr();
    const auto cell_offsets = bins.offsetsPtr();

    const amrex::Real target_ratio = m_target_ratio;

        // Loop over cells
    amrex::ParallelFor( n_cells,
        [=] AMREX_GPU_DEVICE (int i_cell) noexcept
        {
            // The particles that are in the cell `i_cell` are
            // given by the `indices[cell_start:cell_stop]`
            const auto cell_start = cell_offsets[i_cell];
            const auto cell_stop  = cell_offsets[i_cell+1];
            const int cell_numparts = cell_stop - cell_start;

            amrex::Real average_weight = 0._rt;

            for (int i = cell_start; i < cell_stop; ++i)
            {
                average_weight += w[indices[i]]/cell_numparts;
            }

            const amrex::Real level_weight = average_weight*target_ratio;

            amrex::Real random_number;

            for (int i = cell_start; i < cell_stop; ++i)
            {
                if (w[indices[i]] > level_weight) {continue;}

                else
                {
                    random_number = amrex::Random();
                    if (random_number > w[indices[i]]/level_weight)
                    {
                        particle_ptr[indices[i]].id() = -1;
                    }
                    else
                    {
                        w[indices[i]] = level_weight;
                    }
                }
            }
        }
    );
}
