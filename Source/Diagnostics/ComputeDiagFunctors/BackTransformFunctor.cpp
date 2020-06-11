#include "BackTransformFunctor.H"

using namespace amrex;

BackTransformFunctor::BackTransformFunctor(amrex::MultiFab const * mf_src, int lev,
                                           int ncomp, amrex::IntVect crse_ratio)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev)
{}
