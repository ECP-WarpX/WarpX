#include "CellCenterFunctor.H"
#include "Average.H"

using namespace amrex;

CellCenterFunctor::CellCenterFunctor(amrex::MultiFab* mf_src)
    : ComputeDiagFunctor(), m_mf_src(mf_src){}

void
CellCenterFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp, int ncomp) const
{
    Average::ToCellCenter ( mf_dst, *m_mf_src, dcomp, 0, 0, ncomp );
}
