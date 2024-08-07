/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ParticleUtils.H"

#include "WarpX.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Dim3.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

namespace ParticleUtils
{

    using namespace amrex;

    // Define shortcuts for frequently-used type names
    using ParticleType = typename WarpXParticleContainer::ParticleType;
    using ParticleTileType = typename WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = typename ParticleTileType::ParticleTileDataType;
    using ParticleBins = DenseBins<ParticleTileDataType>;
    using index_type = typename ParticleBins::index_type;


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect getParticleCellIndex (const ParticleType& p,
                                  GpuArray<Real,AMREX_SPACEDIM> const& plo,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxi,
                                  Dim3 lo) {

        int ii = 0, jj = 0, kk = 0;
        const ParticleReal x = p.pos(0);
        const Real lx = (x - plo[0]) * dxi[0] - lo.x;
        ii = static_cast<int>(Math::floor(lx));
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
        const ParticleReal y = p.pos(1);
        const Real ly = (y - plo[1]) * dxi[1] - lo.y;
        jj = static_cast<int>(Math::floor(ly));
#endif
#if defined(WARPX_DIM_3D)
        const ParticleReal z = p.pos(2);
        const Real lz = (z - plo[2]) * dxi[2] - lo.z;
        kk = static_cast<int>(Math::floor(lz));
#endif

        return IntVect{AMREX_D_DECL(ii, jj, kk)};
    }

    /* Find the particles and count the particles that are in each cell.
       Note that this does *not* rearrange particle arrays */
    ParticleBins
    findParticlesInEachCell (int lev,
                             MFIter const & mfi,
                             ParticleTileType & ptile) {

        // Extract particle structures for this tile
        int const np = ptile.numParticles();
        auto ptd = ptile.getParticleTileData();

        // Extract box properties
        Geometry const& geom = WarpX::GetInstance().Geom(lev);
        Box const& cbx = mfi.tilebox(IntVect::TheZeroVector()); //Cell-centered box
        const auto lo = lbound(cbx);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        // Find particles that are in each cell;
        // results are stored in the object `bins`.
        ParticleBins bins;
        bins.build(np, ptd, cbx,
            // Pass lambda function that returns the cell index
            [=] AMREX_GPU_DEVICE (ParticleType const & p) noexcept -> IntVect
            {
                return getParticleCellIndex(p, plo, dxi, lo);
            });

        return bins;
    }

} // namespace ParticleUtils
