#include "PartPerGridFunctor.H"
#include "Average.H"

using namespace amrex;

PartPerGridFunctor::PartPerGridFunctor(amrex::MultiFab* mf_src, int lev, int ncomp)
    : ComputeDiagFunctor(ncomp), m_lev(lev)
{
    // mf_src will not be used, let's make sure it's null.
    AMREX_ALWAYS_ASSERT(mf_src == nullptr);
    // Write only in one output component.
    AMREX_ALWAYS_ASSERT(ncomp == 1);
}

void
PartPerGridFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    // Make alias MultiFab* pointing to the component of mf_dst where the
    // number of particles per cell is to be written.
    MultiFab* mf_dst_dcomp = new MultiFab(mf_dst, amrex::make_alias, dcomp, 1);
    // Set value to 0, and increment the value in each cell with ppc.
    mf_dst_dcomp->setVal(0._rt);
    warpx.GetPartContainer().Increment(*mf_dst_dcomp, m_lev);
}
