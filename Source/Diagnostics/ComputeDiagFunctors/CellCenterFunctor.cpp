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
    const amrex::DistributionMapping dm = warpx.DistributionMap(lev);
    const int modes = warpx.n_rz_azimuthal_modes;
    std::unique_ptr<amrex::MultiFab> divE = std::make_unique<amrex::MultiFab>(ba, dm, modes, ng);
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
    const amrex::DistributionMapping dm = warpx.DistributionMap(lev);
    const int modes = warpx.n_rz_azimuthal_modes;
    std::unique_ptr<amrex::MultiFab> divB = std::make_unique<amrex::MultiFab>(ba, dm, modes, ng);
    auto Bfield = warpx.get_array_Bfield_aux(lev);
    warpx.ComputeDivB(*divB, 0, Bfield, WarpX::CellSize(lev));
    return divB.get();
}

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::PartCell> (const int lev, const int dir)
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // Temporary cell-centered, single-component MultiFab for storing particles per cell.
    const amrex::BoxArray& ba = warpx.boxArray(lev);
    const amrex::DistributionMapping dm = warpx.DistributionMap(lev);
    // Set value to 0, and increment the value in each cell with ppc.
    std::unique_ptr<amrex::MultiFab> ppc_mf = std::make_unique<amrex::MultiFab>(ba, dm, 1, ng);
    ppc_mf->setVal(0._rt);
    // Compute ppc which includes a summation over all species.
    warpx.GetPartContainer().Increment(*ppc_mf, lev);
    return ppc_mf.get();
}

template <>
amrex::MultiFab * get_field_pointer <FieldTypes::PartGrid> (const int lev, const int dir)
{
    auto& warpx = WarpX::GetInstance();
    const Vector<long>& npart_in_grid = warpx.GetPartContainer().NumberOfParticlesInGrid(lev);
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // Temporary MultiFab containing number of particles per grid.
    // (stored as constant for all cells in each grid)
    const amrex::BoxArray& ba = warpx.boxArray(lev);
    const amrex::DistributionMapping dm = warpx.DistributionMap(lev);
    std::unique_ptr<amrex::MultiFab> ppg_mf = std::make_unique<amrex::MultiFab>(ba, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(*ppg_mf); mfi.isValid(); ++mfi) {
        (*ppg_mf)[mfi].setVal<RunOn::Host>(static_cast<Real>(npart_in_grid[mfi.index()]));
    }
    return ppg_mf.get();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        const amrex::DistributionMapping dm = warpx.DistributionMap(m_lev);
        MultiFab mf_dst_stag(mf_src->boxArray(), dm, 1, mf_src->nGrowVect());
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
template class CellCenterFunctor<FieldTypes::PartCell>;
template class CellCenterFunctor<FieldTypes::PartGrid>;
