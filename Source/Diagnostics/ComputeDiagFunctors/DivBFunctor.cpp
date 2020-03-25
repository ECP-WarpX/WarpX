#include "DivBFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

DivBFunctor::DivBFunctor(std::array<const amrex::MultiFab*, 3> arr_mf_src, int lev, int ncomp)
    : ComputeDiagFunctor(ncomp), m_arr_mf_src(arr_mf_src), m_lev(lev)
{}

void
DivBFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp) const
{
    WarpX::ComputeDivB(mf_dst, dcomp, m_arr_mf_src, WarpX::CellSize(m_lev) );
}
