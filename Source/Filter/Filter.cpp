/* Copyright 2019-2021 Andrew Myers, Maxence Thevenet, Weiqun Zhang, Axel Huebl
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Filter.H"

#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include <algorithm>

using namespace amrex;

void
Filter::ApplyStencil (MultiFab& dstmf, const MultiFab& srcmf, const int lev, int scomp, int dcomp, int ncomp)
{
    WARPX_PROFILE("Filter::ApplyStencil(MultiFab)");
    ncomp = std::min(ncomp, srcmf.nComp());

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    for (MFIter mfi(dstmf); mfi.isValid(); ++mfi)
    {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        amrex::Real wt = amrex::second();

        const auto& src = srcmf.array(mfi);
        const auto& dst = dstmf.array(mfi);
        const Box& tbx = mfi.growntilebox();

        // Apply filter
        DoFilter(tbx, src, dst, scomp, dcomp, ncomp);

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }
}

void
Filter::ApplyStencil (FArrayBox& dstfab, const FArrayBox& srcfab,
                      const Box& tbx, int scomp, int dcomp, int ncomp)
{
    WARPX_PROFILE("Filter::ApplyStencil(FArrayBox)");
    ncomp = std::min(ncomp, srcfab.nComp());
    const auto& src = srcfab.array();
    const auto& dst = dstfab.array();

    // Apply filter
    DoFilter(tbx, src, dst, scomp, dcomp, ncomp);
}

void Filter::DoFilter (const Box& tbx,
                       Array4<Real const> const& src,
                       Array4<Real      > const& dst,
                       int scomp, int dcomp, int ncomp)
{
    amrex::Real const* AMREX_RESTRICT sx = stencil_x.data();
#if (AMREX_SPACEDIM == 3)
    amrex::Real const* AMREX_RESTRICT sy = stencil_y.data();
#endif
    amrex::Real const* AMREX_RESTRICT sz = stencil_z.data();
    Dim3 slen_local = slen;

#if (AMREX_SPACEDIM == 3)
    AMREX_PARALLEL_FOR_4D( tbx, ncomp, i, j, k, n,
    {
        Real d = 0.0;

        // Pad source array with zeros beyond ghost cells
        // for out-of-bound accesses due to large-stencil operations
        const auto src_zeropad = [src] (const int jj, const int kk, const int ll, const int nn) noexcept
        {
            return src.contains(jj,kk,ll) ? src(jj,kk,ll,nn) : 0.0_rt;
        };

        for         (int iz=0; iz < slen_local.z; ++iz){
            for     (int iy=0; iy < slen_local.y; ++iy){
                for (int ix=0; ix < slen_local.x; ++ix){
                    Real sss = sx[ix]*sy[iy]*sz[iz];
                    d += sss*( src_zeropad(i-ix,j-iy,k-iz,scomp+n)
                              +src_zeropad(i+ix,j-iy,k-iz,scomp+n)
                              +src_zeropad(i-ix,j+iy,k-iz,scomp+n)
                              +src_zeropad(i+ix,j+iy,k-iz,scomp+n)
                              +src_zeropad(i-ix,j-iy,k+iz,scomp+n)
                              +src_zeropad(i+ix,j-iy,k+iz,scomp+n)
                              +src_zeropad(i-ix,j+iy,k+iz,scomp+n)
                              +src_zeropad(i+ix,j+iy,k+iz,scomp+n));
                }
            }
        }

        dst(i,j,k,dcomp+n) = d;
    });
#else
    AMREX_PARALLEL_FOR_4D( tbx, ncomp, i, j, k, n,
    {
        Real d = 0.0;

        // Pad source array with zeros beyond ghost cells
        // for out-of-bound accesses due to large-stencil operations
        const auto src_zeropad = [src] (const int jj, const int kk, const int ll, const int nn) noexcept
        {
            return src.contains(jj,kk,ll) ? src(jj,kk,ll,nn) : 0.0_rt;
        };

        for         (int iz=0; iz < slen_local.z; ++iz){
            for     (int iy=0; iy < slen_local.y; ++iy){
                for (int ix=0; ix < slen_local.x; ++ix){
                    Real sss = sx[ix]*sz[iy];
                    d += sss*( src_zeropad(i-ix,j-iy,k,scomp+n)
                              +src_zeropad(i+ix,j-iy,k,scomp+n)
                              +src_zeropad(i-ix,j+iy,k,scomp+n)
                              +src_zeropad(i+ix,j+iy,k,scomp+n));
                }
            }
        }

        dst(i,j,k,dcomp+n) = d;
    });
#endif
}
