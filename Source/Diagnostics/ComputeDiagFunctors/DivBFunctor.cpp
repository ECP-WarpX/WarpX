#include "DivBFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

DivBFunctor::DivBFunctor(const std::array<const amrex::MultiFab* const, 3> arr_mf_src, const int lev, amrex::IntVect diag_crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, diag_crse_ratio), m_arr_mf_src(arr_mf_src), m_lev(lev)
{}

void
DivBFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp) const
{
    WarpX::ComputeDivB(mf_dst, dcomp, m_arr_mf_src, WarpX::CellSize(m_lev) );
}
