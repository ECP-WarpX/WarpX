/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXCommUtil.H"

#include <AMReX.H>
#include <AMReX_BaseFab.H>
#include <AMReX_IntVect.H>
#include <AMReX_FabArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

namespace WarpXCommUtil {

void ParallelCopy (amrex::MultiFab&            dst,
                   const amrex::MultiFab&      src,
                   int                         src_comp,
                   int                         dst_comp,
                   int                         num_comp,
                   const amrex::IntVect&       src_nghost,
                   const amrex::IntVect&       dst_nghost,
                   const amrex::Periodicity&   period,
                   amrex::FabArrayBase::CpOp   op)
{
    BL_PROFILE("WarpXCommUtil::ParallelCopy");

    using WarpXCommUtil::comm_float_type;

    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > src_tmp(src.boxArray(),
                                                                  src.DistributionMap(),
                                                                  num_comp,
                                                                  src_nghost);
        mixedCopy(src_tmp, src, src_comp, 0, num_comp, src_nghost);

        amrex::FabArray<amrex::BaseFab<comm_float_type> > dst_tmp(dst.boxArray(),
                                                                  dst.DistributionMap(),
                                                                  num_comp,
                                                                  dst_nghost);

        mixedCopy(dst_tmp, dst, dst_comp, 0, num_comp, dst_nghost);

        dst_tmp.ParallelCopy(src_tmp, 0, 0, num_comp,
                             src_nghost, dst_nghost, period, op);

        mixedCopy(dst, dst_tmp, 0, dst_comp, num_comp, dst_nghost);
    }
    else
    {
        dst.ParallelCopy(src, src_comp, dst_comp, num_comp, src_nghost, dst_nghost, period, op);
    }
}

void ParallelAdd (amrex::MultiFab&            dst,
                  const amrex::MultiFab&      src,
                  int                         src_comp,
                  int                         dst_comp,
                  int                         num_comp,
                  const amrex::IntVect&       src_nghost,
                  const amrex::IntVect&       dst_nghost,
                  const amrex::Periodicity&   period)
{
    WarpXCommUtil::ParallelCopy(dst, src, src_comp, dst_comp, num_comp, src_nghost, dst_nghost, period,
                                amrex::FabArrayBase::ADD);
}

void FillBoundary (amrex::MultiFab& mf, const amrex::Periodicity& period)
{
    BL_PROFILE("WarpXCommUtil::FillBoundary");

    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > mf_tmp(mf.boxArray(),
                                                                 mf.DistributionMap(),
                                                                 mf.nComp(),
                                                                 mf.nGrowVect());

        mixedCopy(mf_tmp, mf, 0, 0, mf.nComp(), mf.nGrowVect());

        mf_tmp.FillBoundary(period);

        mixedCopy(mf, mf_tmp, 0, 0, mf.nComp(), mf.nGrowVect());
    }
    else
    {
        mf.FillBoundary(period);
    }
}

void FillBoundary (amrex::MultiFab&          mf,
                   amrex::IntVect            ng,
                   const amrex::Periodicity& period)
{
    BL_PROFILE("WarpXCommUtil::FillBoundary");

    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > mf_tmp(mf.boxArray(),
                                                            mf.DistributionMap(),
                                                            mf.nComp(),
                                                            mf.nGrowVect());

        mixedCopy(mf_tmp, mf, 0, 0, mf.nComp(), mf.nGrowVect());

        mf_tmp.FillBoundary(ng, period);

        mixedCopy(mf, mf_tmp, 0, 0, mf.nComp(), mf.nGrowVect());
    }
    else
    {
        mf.FillBoundary(ng, period);
    }
}

void FillBoundary (amrex::iMultiFab& imf, const amrex::Periodicity& period)
{
    BL_PROFILE("WarpXCommUtil::FillBoundary");

    imf.FillBoundary(period);
}

void FillBoundary (amrex::iMultiFab&         imf,
                   amrex::IntVect            ng,
                   const amrex::Periodicity& period)
{
    BL_PROFILE("WarpXCommUtil::FillBoundary");
    imf.FillBoundary(ng, period);
}

void
FillBoundary (amrex::Vector<amrex::MultiFab*> const& mf, const amrex::Periodicity& period)
{
    for (auto x : mf) {
        WarpXCommUtil::FillBoundary(*x, period);
    }
}

void SumBoundary (amrex::MultiFab& mf, const amrex::Periodicity& period)
{
    BL_PROFILE("WarpXCommUtil::SumBoundary");

    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > mf_tmp(mf.boxArray(),
                                                                 mf.DistributionMap(),
                                                                 mf.nComp(),
                                                                 mf.nGrowVect());

        mixedCopy(mf_tmp, mf, 0, 0, mf.nComp(), mf.nGrowVect());

        mf_tmp.SumBoundary(period);

        mixedCopy(mf, mf_tmp, 0, 0, mf.nComp(), mf.nGrowVect());
    }
    else
    {
        mf.SumBoundary(period);
    }
}

void SumBoundary (amrex::MultiFab&          mf,
                  int                       start_comp,
                  int                       num_comps,
                  amrex::IntVect            ng,
                  const amrex::Periodicity& period)
{
    BL_PROFILE("WarpXCommUtil::SumBoundary");

    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > mf_tmp(mf.boxArray(),
                                                                 mf.DistributionMap(),
                                                                 num_comps,
                                                                 ng);
        mixedCopy(mf_tmp, mf, start_comp, 0, num_comps, ng);

        mf_tmp.SumBoundary(0, num_comps, ng, period);

        mixedCopy(mf, mf_tmp, 0, start_comp, num_comps, ng);
    }
    else
    {
        mf.SumBoundary(start_comp, num_comps, ng, period);
    }
}

void SumBoundary (amrex::MultiFab&          mf,
                  int                       start_comp,
                  int                       num_comps,
                  amrex::IntVect            src_ng,
                  amrex::IntVect            dst_ng,
                  const amrex::Periodicity& period)
{
    BL_PROFILE("WarpXCommUtil::SumBoundary");

    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > mf_tmp(mf.boxArray(),
                                                                 mf.DistributionMap(),
                                                                 num_comps,
                                                                 mf.nGrowVect());
        mixedCopy(mf_tmp, mf, start_comp, 0, num_comps, mf.nGrowVect());

        mf_tmp.SumBoundary(0, num_comps, src_ng, dst_ng, period);

        mixedCopy(mf, mf_tmp, 0, start_comp, num_comps, dst_ng);
    }
    else
    {
        mf.SumBoundary(start_comp, num_comps, src_ng, dst_ng, period);
    }
}

void OverrideSync (amrex::MultiFab&          mf,
                   const amrex::Periodicity& period)
{
    if (mf.ixType().cellCentered()) return;

    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > mf_tmp(mf.boxArray(),
                                                                 mf.DistributionMap(),
                                                                 mf.nComp(),
                                                                 mf.nGrowVect());

        mixedCopy(mf_tmp, mf, 0, 0, mf.nComp(), mf.nGrowVect());

        auto msk = mf.OwnerMask(period);
        amrex::OverrideSync(mf_tmp, *msk, period);

        mixedCopy(mf, mf_tmp, 0, 0, mf.nComp(), mf.nGrowVect());
    }
    else
    {
        mf.OverrideSync(period);
    }
}

}
