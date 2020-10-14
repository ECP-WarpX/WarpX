#include "WarpXCommUtil.H"

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
    if (WarpX::do_single_precision_comms)
    {
        amrex::FabArray<amrex::BaseFab<comm_float_type> > src_tmp(src.boxArray(),
                                                                  src.DistributionMap(),
                                                                  src.nComp(),
                                                                  src.nGrowVect());
        mixedCopy(src_tmp, src, 0, 0, src.nComp(), src.nGrowVect());

        amrex::FabArray<amrex::BaseFab<comm_float_type> > dst_tmp(dst.boxArray(),
                                                                  dst.DistributionMap(),
                                                                  dst.nComp(),
                                                                  dst.nGrowVect());

        mixedCopy(dst_tmp, dst, 0, 0, dst.nComp(), dst.nGrowVect());

        dst_tmp.ParallelCopy(src_tmp, src_comp, dst_comp, num_comp,
                             src_nghost, dst_nghost, period, op);

        mixedCopy(dst, dst_tmp, 0, 0, dst.nComp(), dst.nGrowVect());
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
    ParallelCopy(dst, src, src_comp, dst_comp, num_comp, src_nghost, dst_nghost, period,
                 amrex::FabArrayBase::ADD);
}

void FillBoundary (amrex::MultiFab& mf, const amrex::Periodicity& period)
{
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

void
FillBoundary (amrex::Vector<amrex::MultiFab*> const& mf, const amrex::Periodicity& period)
{
    for (auto x : mf) {
        WarpXCommUtil::FillBoundary(*x, period);
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
