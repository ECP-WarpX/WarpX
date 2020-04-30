#include "CellCenterFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::E>::get_field_diag (const int lev, const int dir) const
{
    auto & warpx = WarpX::GetInstance();
    m_do_delete = false;
    return warpx.get_pointer_Efield_aux(lev, dir);
};

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::B>::get_field_diag (const int lev, const int dir) const
{
    auto & warpx = WarpX::GetInstance();
    m_do_delete = false;
    return warpx.get_pointer_Bfield_aux(lev, dir);
};

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::J>::get_field_diag (const int lev, const int dir) const
{
    auto & warpx = WarpX::GetInstance();
    m_do_delete = false;
    return warpx.get_pointer_current_fp(lev, dir);
};

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::rho>::get_field_diag (const int lev, const int dir) const
{
    auto & warpx = WarpX::GetInstance();
#ifdef WARPX_USE_PSATD
    // rho_new is stored in component 1 of rho_fp when using PSATD
    amrex::MultiFab * rho_new = new amrex::MultiFab(*warpx.get_pointer_rho_fp(lev), amrex::make_alias, 1, 1);
    m_do_delete = true;
    return rho_new;
#else
    m_do_delete = false;
    return warpx.get_pointer_rho_fp(lev);
#endif
};

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::F>::get_field_diag (const int lev, const int dir) const
{
    auto & warpx = WarpX::GetInstance();
    m_do_delete = false;
    return warpx.get_pointer_F_fp(lev);
};


template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::divE>::get_field_diag (const int lev, const int dir) const
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
    amrex::MultiFab * divE = new amrex::MultiFab(ba, dm, modes, ng);
    warpx.ComputeDivE(*divE, lev);
    m_do_delete = true;
    return divE;
}

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::divB>::get_field_diag (const int lev, const int dir) const
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
    amrex::MultiFab * divB = new amrex::MultiFab(ba, dm, modes, ng);
    auto Bfield = warpx.get_array_Bfield_aux(lev);
    warpx.ComputeDivB(*divB, 0, Bfield, WarpX::CellSize(lev));
    m_do_delete = true;
    return divB;
}

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::PartCell>::get_field_diag (const int lev, const int dir) const
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
    amrex::MultiFab * ppc_mf = new amrex::MultiFab(ba, dm, 1, ng);
    ppc_mf->setVal(0._rt);
    // Compute ppc which includes a summation over all species.
    warpx.GetPartContainer().Increment(*ppc_mf, lev);
    m_do_delete = true;
    return ppc_mf;
}

template <>
amrex::MultiFab *
CellCenterFunctor<FieldTypes::PartGrid>::get_field_diag (const int lev, const int dir) const
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
    amrex::MultiFab * ppg_mf = new amrex::MultiFab(ba, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(*ppg_mf); mfi.isValid(); ++mfi) {
        (*ppg_mf)[mfi].setVal<RunOn::Host>(static_cast<Real>(npart_in_grid[mfi.index()]));
    }
    m_do_delete = true;
    return ppg_mf;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <FieldTypes::TypeEnum FIELDTYPE>
void
CellCenterFunctor<FIELDTYPE>::operator()(amrex::MultiFab& mf_dst, int dcomp, const amrex::IntVect crse_ratio) const
{

    MultiFab const * const mf_src = get_field_diag(m_lev, m_dir);

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
    if (m_do_delete) {
        delete mf_src;
    }
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
