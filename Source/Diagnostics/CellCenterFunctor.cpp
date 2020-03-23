#include "CellCenterFunctor.H"
#include "Average.H"

using namespace amrex;

CellCenterFunctor::CellCenterFunctor(amrex::MultiFab* mf_src, int lev, bool mode_avg, int ncomp)
    : ComputeDiagFunctor(ncomp), m_mf_src(mf_src), m_lev(lev), m_mode_avg(mode_avg)
{}

void
CellCenterFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp) const
{
#ifdef WARPX_DIM_RZ
    if (m_mode_avg) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into 1 single component");
        // In cylindrical gemoetry, sum real part of all modes in m_mf_src in
        // temporary multifab mf_dst_stag, and cell-center it to mf_dst.
        auto& warpx = WarpX::GetInstance();
        MultiFab mf_dst_stag(m_mf_src->boxArray(), warpx.DistributionMap(m_lev), 1, m_mf_src->nGrowVect());
        // Mode 0
        MultiFab::Copy(mf_dst_stag, *m_mf_src, 0, 0, 1, m_mf_src->nGrowVect());
        for (int ic=1 ; ic < m_mf_src->nComp() ; ic += 2) {
            // All modes > 0
            MultiFab::Add(mf_dst_stag, *m_mf_src, ic, 0, 1, m_mf_src->nGrowVect());
        }
        Average::ToCellCenter ( mf_dst, mf_dst_stag, dcomp, 0, 0, nComp() );
    } else {
        Average::ToCellCenter ( mf_dst, *m_mf_src, dcomp, 0, 0, nComp() );
    }
#else
    // In cartesian geometry, cell-center m_mf_src to mf_dst.
    Average::ToCellCenter ( mf_dst, *m_mf_src, dcomp, 0, 0, nComp() );
#endif
}
