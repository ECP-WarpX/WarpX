#include "CellCenterFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

template <FieldTypes::TypeEnum FIELDTYPE>
amrex::MultiFab * get_field_pointer (const int lev, const int dir)
{
    return nullptr;
};

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::E> (const int lev, const int dir)
{
    auto & warpx = WarpX::GetInstance();
    return warpx.get_pointer_Efield_aux(lev, dir);
};

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::B> (const int lev, const int dir)
{
    auto & warpx = WarpX::GetInstance();
    return warpx.get_pointer_Bfield_aux(lev, dir);
};

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::J> (const int lev, const int dir)
{
    auto & warpx = WarpX::GetInstance();
    return warpx.get_pointer_current_fp(lev, dir);
};

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::rho> (const int lev, const int dir)
{
    auto & warpx = WarpX::GetInstance();
#ifdef WARPX_USE_PSATD
    // rho_new is stored in component 1 of rho_fp when using PSATD
    std::unique_ptr<amrex::MultiFab> rho_new = std::make_unique<amrex::MultiFab>(*warpx.get_pointer_rho_fp(lev), amrex::make_alias, 1, 1);
    return rho_new.get();
#else
    return warpx.get_pointer_rho_fp(lev);
#endif
};

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::F> (const int lev, const int dir)
{
    auto & warpx = WarpX::GetInstance();
    return warpx.get_pointer_F_fp(lev);
};


template <>
amrex::MultiFab * get_field_pointer <FieldTypes::divE> (const int lev, const int dir)
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // For staggered and nodal calculations, divE is computed on the nodes.
    // The temporary divE MultiFab is generated to comply with the location of divE.
    const amrex::BoxArray& ba = amrex::convert(warpx.boxArray(lev),amrex::IntVect::TheUnitVector());
    std::unique_ptr<amrex::MultiFab> divE = std::make_unique<amrex::MultiFab>(ba, warpx.DistributionMap(lev), 1, ng);
    warpx.ComputeDivE(*divE, lev);
    return divE.get();
}

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::divB> (const int lev, const int dir)
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // A cell-centered divB multifab spanning the entire domain is generated
    // and divB is computed on the cell-center, with ng=1.
    const amrex::BoxArray& ba = warpx.boxArray(lev);
    std::unique_ptr<amrex::MultiFab> divB = std::make_unique<amrex::MultiFab>(ba, warpx.DistributionMap(lev), 1, ng);
    auto Bfield = warpx.get_array_Bfield_aux(lev);
    warpx.ComputeDivB(*divB, 0, Bfield, WarpX::CellSize(lev));
    return divB.get();
}

template <FieldTypes::TypeEnum FIELDTYPE>
CellCenterFunctor<FIELDTYPE>::CellCenterFunctor (const int lev, const int dir,
                                                 const bool convertRZmodes2cartesian, const int ncomp)
    : ComputeDiagFunctor(ncomp), m_lev(lev), m_dir(dir),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{}

template <FieldTypes::TypeEnum FIELDTYPE>
void
CellCenterFunctor<FIELDTYPE>::operator()(amrex::MultiFab& mf_dst, int dcomp, const amrex::IntVect crse_ratio) const
{

    MultiFab const * const mf_src = get_field_pointer <FIELDTYPE> (m_lev, m_dir);

#ifdef WARPX_DIM_RZ
    if (m_convertRZmodes2cartesian) {
        // In cylindrical geometry, sum real part of all modes of mf_src in
        // temporary multifab mf_dst_stag, and cell-center it to mf_dst.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into 1 single component");
        auto& warpx = WarpX::GetInstance();
        MultiFab mf_dst_stag(mf_src->boxArray(), warpx.DistributionMap(m_lev), 1, mf_src->nGrowVect());
        // Mode 0
        MultiFab::Copy(mf_dst_stag, *mf_src, 0, 0, 1, mf_src->nGrowVect());
        for (int ic=1 ; ic < mf_src->nComp() ; ic += 2) {
            // All modes > 0
            MultiFab::Add(mf_dst_stag, *mf_src, ic, 0, 1, mf_src->nGrowVect());
        }
        Average::CoarsenAndInterpolate(mf_dst, mf_dst_stag, dcomp, 0, nComp(), 0,  crse_ratio);
    } else {
        Average::CoarsenAndInterpolate(mf_dst, *mf_src, dcomp, 0, nComp(), 0, crse_ratio);
    }
#else
    // In cartesian geometry, coarsen and interpolate from simulation MultiFab, mf_src,
    // to output diagnostic MultiFab, mf_dst.
    Average::CoarsenAndInterpolate(mf_dst, *mf_src, dcomp, 0, nComp(), 0, crse_ratio);
#endif
}

// These are needed to force the compiler to build these versions
template class CellCenterFunctor<FieldTypes::E>;
template class CellCenterFunctor<FieldTypes::B>;
template class CellCenterFunctor<FieldTypes::J>;
template class CellCenterFunctor<FieldTypes::rho>;
template class CellCenterFunctor<FieldTypes::F>;
template class CellCenterFunctor<FieldTypes::divE>;
template class CellCenterFunctor<FieldTypes::divB>;
