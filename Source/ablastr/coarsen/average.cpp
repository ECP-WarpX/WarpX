/* Copyright 2022 Edoardo Zoni, Remi Lehe, Prabhat Kumar
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "average.H"

#include "ablastr/utils/TextMsg.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace ablastr::coarsen::average
{
    void
    Loop (
        amrex::MultiFab & mf_dst,
        amrex::MultiFab const & mf_src,
        int const ncomp,
        amrex::IntVect const ngrow,
        amrex::IntVect const crse_ratio
    )
    {
        // Staggering of source fine MultiFab and destination coarse MultiFab
        amrex::IntVect const stag_src = mf_src.boxArray().ixType().toIntVect();
        amrex::IntVect const stag_dst = mf_dst.boxArray().ixType().toIntVect();

        // Auxiliary integer arrays (always 3D)
        auto sf = amrex::GpuArray<int,3>{0,0,0}; // staggering of source fine MultiFab
        auto sc = amrex::GpuArray<int,3>{0,0,0}; // staggering of destination coarse MultiFab
        auto cr = amrex::GpuArray<int,3>{1,1,1}; // coarsening ratio
        for (int i=0; i<AMREX_SPACEDIM; ++i)
        {
             sf[i] = stag_src[i];
             sc[i] = stag_dst[i];
             cr[i] = crse_ratio[i];
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        // Loop over boxes (or tiles if not on GPU)
        for (amrex::MFIter mfi(mf_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Tiles defined at the coarse level
            amrex::Box const & bx = mfi.growntilebox(ngrow);
            amrex::Array4<amrex::Real> const &arr_dst = mf_dst.array(mfi);
            amrex::Array4<amrex::Real const> const &arr_src = mf_src.const_array(mfi);
            ParallelFor(bx, ncomp,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                            arr_dst(i, j, k, n) = Interp(
                                arr_src, sf, sc, cr, i, j, k, n);
                        });
        }
    }

    void
    Coarsen (
        amrex::MultiFab & mf_dst,
        amrex::MultiFab const & mf_src,
        amrex::IntVect const crse_ratio
    )
    {
        BL_PROFILE("ablastr::coarsen::Coarsen()");

        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(mf_src.ixType() == mf_dst.ixType(),
                                           "source MultiFab and destination MultiFab have different IndexType");

        // Number of guard cells to fill on coarse patch and number of components
        const amrex::IntVect ngrow = (mf_src.nGrowVect() + crse_ratio-1) / crse_ratio; // round up int division
        const int ncomp = mf_src.nComp();

        Loop(mf_dst, mf_src, ncomp, ngrow, crse_ratio);
    }

    void
    CoarseningPointsAndWeights (
        amrex::Vector<amrex::Real> &weights_x,
        amrex::Vector<amrex::Real> &weights_y,
        amrex::Vector<amrex::Real> &weights_z,
        amrex::GpuArray<int, 3> &src_index_min,
        amrex::GpuArray<int, 3> const &coarsen_ratio,
        amrex::GpuArray<int, 3> const &stag_fine_src,
        amrex::GpuArray<int, 3> const &stag_crse_des
    )
    {
        using namespace amrex::literals;

        int num_points[3];
        bool useHalf[3];

        for (int l = 0; l < 3; ++l) {
            int const twoImin = -coarsen_ratio[l]*stag_crse_des[l] + stag_fine_src[l] - 1;
            if (twoImin % 2 == 0) {
                src_index_min[l] = twoImin/2;
                num_points[l] = coarsen_ratio[l]+1;
                useHalf[l] = true;
            } else {
                src_index_min[l] = (twoImin+1)/2;
                num_points[l] = coarsen_ratio[l];
                useHalf[l] = false;
            }
        }

        // amrex::Real const coarsen_factor = 1.0_rt / static_cast<amrex::Real>(coarsen_ratio[0]*coarsen_ratio[1]*coarsen_ratio[2]);

        weights_x.resize(num_points[0], 1.0_rt);
        weights_y.resize(num_points[1], 1.0_rt);
        weights_z.resize(num_points[2], 1.0_rt);

        if (useHalf[0]) {
            weights_x[0] *= 0.5_rt;
            weights_x[num_points[0]-1] *= 0.5_rt;
        }
        if (useHalf[1]) {
            weights_y[0] *= 0.5_rt;
            weights_y[num_points[1]-1] *= 0.5_rt;
        }
        if (useHalf[2]) {
            weights_z[0] *= 0.5_rt;
            weights_z[num_points[2]-1] *= 0.5_rt;
        }

    }

} // namespace ablastr::coarsen::average
