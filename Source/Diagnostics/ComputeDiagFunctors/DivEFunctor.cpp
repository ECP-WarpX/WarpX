#include "DivEFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

DivEFunctor::DivEFunctor(const std::array<const amrex::MultiFab* const, 3> arr_mf_src, const int lev, const amrex::IntVect diag_crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, diag_crse_ratio), m_arr_mf_src(arr_mf_src), m_lev(lev)
{}

void
DivEFunctor::operator()(amrex::MultiFab& mf_dst, const int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    const int ng = 1;
    // For staggered and nodal calculations, divE is computed on the nodes.
    // The temporary divE multifab is generated to comply with the location of divE
    const BoxArray& ba = amrex::convert(warpx.boxArray(m_lev),IntVect::TheUnitVector());
    MultiFab divE(ba, warpx.DistributionMap(m_lev), 1, ng );
    warpx.ComputeDivE(divE, m_lev);
    // Coarsen and interpolate from divE on the entire domain to cell-centered mf_dst
    Average::CoarsenAndInterpolate(mf_dst, divE, dcomp, 0, nComp(), m_crse_ratio);

}
