#include "PartPerCellFunctor.H"
#include "WarpX.H"
#include "Utils/Average.H"

using namespace amrex;

PartPerCellFunctor::PartPerCellFunctor(const amrex::MultiFab* mf_src, const int lev, amrex::IntVect diag_crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, diag_crse_ratio), m_lev(lev)
{
    // mf_src will not be used, let's make sure it's null.
    AMREX_ALWAYS_ASSERT(mf_src == nullptr);
    // Write only in one output component.
    AMREX_ALWAYS_ASSERT(ncomp == 1);
}

void
PartPerCellFunctor::operator()(amrex::MultiFab& mf_dst, const int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    const int ng = 1;
    // Create a tmp cell-centered, single-component MultiFab to compute ppc
    MultiFab ppc_mf(warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 1, ng);
    // Set value to 0, and increment the value in each cell with ppc.
    ppc_mf.setVal(0._rt);
    // Compute ppc which includes a summation over all species
    warpx.GetPartContainer().Increment(ppc_mf, m_lev);
    // Coarsen and interpolate from ppc_mf to mf_dst
    Average::CoarsenAndInterpolate(mf_dst, ppc_mf, dcomp, 0, nComp(), 0, m_crse_ratio);
}
