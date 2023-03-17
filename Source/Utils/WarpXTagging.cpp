/* Copyright 2019 Axel Huebl, Maxence Thevenet, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <WarpX.H>

#include <AMReX_BaseFab.H>
#include <AMReX_Config.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_SPACE.H>
#include <AMReX_TagBox.H>

#include <AMReX_BaseFwd.H>

using namespace amrex;

void
WarpX::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx = Geom(lev).CellSizeArray();

    const auto ftlo = fine_tag_lo;
    const auto fthi = fine_tag_hi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tags); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.fabbox();
        const auto& fab = tags.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            RealVect pos {AMREX_D_DECL((i+0.5_rt)*dx[0]+problo[0],
                                       (j+0.5_rt)*dx[1]+problo[1],
                                       (k+0.5_rt)*dx[2]+problo[2])};
            if (pos > ftlo && pos < fthi) {
                fab(i,j,k) = TagBox::SET;
            }
        });
    }
}
