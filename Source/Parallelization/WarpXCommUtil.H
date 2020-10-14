#ifndef WARPX_COMMUTIL_H_
#define WARPX_COMMUTIL_H_

#include "AMReX_FabArray.H"
#include "AMReX_Periodicity.H"
#include "AMReX_TypeTraits.H"

#include "WarpX.H"

namespace WarpXCommUtil
{

    using comm_float_type = float;

template <class FAB1, class FAB2>
void
mixedCopy (amrex::FabArray<FAB1>& dst, amrex::FabArray<FAB2> const& src, int srccomp, int dstcomp, int numcomp, const amrex::IntVect& nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(dst,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto const srcFab = src.array(mfi);
            auto       dstFab = dst.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FUSIBLE ( bx, numcomp, i, j, k, n,
            {
                dstFab(i,j,k,dstcomp+n) = srcFab(i,j,k,srccomp+n);
            });
        }
    }
}

void ParallelCopy (amrex::MultiFab&            dst,
                   const amrex::MultiFab&      src,
                   int                         src_comp,
                   int                         dst_comp,
                   int                         num_comp,
                   const amrex::IntVect&       src_nghost,
                   const amrex::IntVect&       dst_nghost,
                   const amrex::Periodicity&   period = amrex::Periodicity::NonPeriodic(),
                   amrex::FabArrayBase::CpOp   op = amrex::FabArrayBase::COPY);

void ParallelAdd (amrex::MultiFab&            dst,
                  const amrex::MultiFab&      src,
                  int                         src_comp,
                  int                         dst_comp,
                  int                         num_comp,
                  const amrex::IntVect&       src_nghost,
                  const amrex::IntVect&       dst_nghost,
                  const amrex::Periodicity&   period = amrex::Periodicity::NonPeriodic());

void FillBoundary (amrex::MultiFab&          mf,
                   const amrex::Periodicity& period = amrex::Periodicity::NonPeriodic());

void FillBoundary (amrex::MultiFab&          mf,
                   amrex::IntVect            ng,
                   const amrex::Periodicity& period = amrex::Periodicity::NonPeriodic());
void
FillBoundary (amrex::Vector<amrex::MultiFab*> const& mf, const amrex::Periodicity& period);

void OverrideSync (amrex::MultiFab&          mf,
                   const amrex::Periodicity& period = amrex::Periodicity::NonPeriodic());
}

#endif
