/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldSolver/ImplicitSolvers/WarpXSolverVec.H"

void WarpXSolverVec::SetDotMask( const amrex::Vector<amrex::Geometry>&  a_Geom )
{
    if (m_dot_mask_defined) { return; }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            IsDefined(),
        "WarpXSolverVec::SetDotMask() called from undefined instance ");

    m_dotMask.resize(m_num_amr_levels);
    for ( int n = 0; n < 3; n++) {
        const amrex::BoxArray& grids = m_field_vec[0][n]->boxArray();
        const amrex::MultiFab tmp( grids, m_field_vec[0][n]->DistributionMap(),
                                   1, 0, amrex::MFInfo().SetAlloc(false) );
        const amrex::Periodicity& period = a_Geom[0].periodicity();
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
        a_X.IsDefined(),
        "WarpXSolverVec::dotProduct(a_X) called with undefined a_X");
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
