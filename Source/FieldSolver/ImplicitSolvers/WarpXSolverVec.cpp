/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldSolver/ImplicitSolvers/WarpXSolverVec.H"
#include "WarpX.H"

void WarpXSolverVec::Define ( WarpX*                    a_WarpX,
                        const warpx::fields::FieldType  a_field_type )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !IsDefined(),
        "WarpXSolverVec::Define(a_vec, a_type) called on already defined WarpXSolverVec");

    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& this_vec(
                                                            a_WarpX->getMultiLevelField(a_field_type));
    m_field_vec.resize(m_num_amr_levels);
    for (int lev=0; lev<m_num_amr_levels; ++lev) {
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& mf_model = *this_vec[lev][n];
            m_field_vec[lev][n] = std::make_unique<amrex::MultiFab>( mf_model.boxArray(),
                                                                     mf_model.DistributionMap(),
                                                                     mf_model.nComp(),
                                                                     amrex::IntVect::TheZeroVector() );
        }
    }
    m_field_type = a_field_type;
    m_is_defined = true;
    SetWarpXPointer(a_WarpX);
    SetDotMask();
}

void WarpXSolverVec::Copy ( const warpx::fields::FieldType  a_field_type )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(),
        "WarpXSolverVec::Copy() called on undefined WarpXSolverVec");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_field_type==m_field_type,
        "WarpXSolverVec::Copy() called with vecs of different types");

    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& this_vec(
                                               m_WarpX->getMultiLevelField(a_field_type));
    for (int lev=0; lev<m_num_amr_levels; ++lev) {
        for (int n=0; n<3; ++n) {
            amrex::MultiFab::Copy( *m_field_vec[lev][n], *this_vec[lev][n], 0, 0, m_ncomp,
                                   amrex::IntVect::TheZeroVector() );
        }
    }
}

void WarpXSolverVec::SetWarpXPointer( WarpX*  a_WarpX )
{
    if (m_warpx_ptr_defined) { return; }
    m_WarpX = a_WarpX;
    m_warpx_ptr_defined = true;
}

void WarpXSolverVec::SetDotMask()
{
    if (m_dot_mask_defined) { return; }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(),
        "WarpXSolverVec::SetDotMask() called from undefined instance ");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_warpx_ptr_defined,
        "WarpXSolverVec::SetDotMask() called before SetWarpXPointer() ");

    const amrex::Vector<amrex::Geometry>& Geom = m_WarpX->Geom();
    m_dotMask.resize(m_num_amr_levels);
    for ( int n = 0; n < 3; n++) {
        const amrex::BoxArray& grids = m_field_vec[0][n]->boxArray();
        const amrex::MultiFab tmp( grids, m_field_vec[0][n]->DistributionMap(),
                                   1, 0, amrex::MFInfo().SetAlloc(false) );
        const amrex::Periodicity& period = Geom[0].periodicity();
        m_dotMask[0][n] = tmp.OwnerMask(period);
    }
    m_dot_mask_defined = true;

    // If the function below is not called, then the following
    // error message occurs after the simulation finishes:
    // malloc_consolidate(): unaligned fastbin chunk detected
    amrex::ExecOnFinalize(WarpXSolverVec::clearDotMask);
}

[[nodiscard]] amrex::Real WarpXSolverVec::dotProduct ( const WarpXSolverVec&  a_X ) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_dot_mask_defined,
        "WarpXSolverVec::dotProduct called with m_dotMask not yet defined");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_X.m_field_type==m_field_type,
        "WarpXSolverVec::dotProduct(X) called with solver vecs of different types");

    amrex::Real result = 0.0;
    const int lev = 0;
    const bool local = true;
    for (int n = 0; n < 3; ++n) {
        auto rtmp = amrex::MultiFab::Dot( *m_dotMask[lev][n],
                                          *m_field_vec[lev][n], 0,
                                          *a_X.getVec()[lev][n], 0, 1, 0, local);
        result += rtmp;
    }
    amrex::ParallelAllReduce::Sum(result, amrex::ParallelContext::CommunicatorSub());
    return result;
}
