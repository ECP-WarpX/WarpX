/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldSolver/ImplicitSolvers/WarpXSolverVec.H"
#include "WarpX.H"

using namespace warpx::fields;

void WarpXSolverVec::Define ( WarpX*     a_WarpX,
                              FieldType  a_array_type,
                              FieldType  a_scalar_type )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !IsDefined(),
        "WarpXSolverVec::Define(a_vec, a_type) called on already defined WarpXSolverVec");

    m_array_type = a_array_type;
    m_scalar_type = a_scalar_type;

    m_array_vec.resize(m_num_amr_levels);
    m_scalar_vec.resize(m_num_amr_levels);

    // Define the 3D vector field data container
    if (m_array_type != FieldType::None) {

        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n = 0; n < 3; n++) {
                const amrex::MultiFab* mf_model = a_WarpX->getFieldPointer(m_array_type,lev,n);
                m_array_vec[lev][n] = std::make_unique<amrex::MultiFab>( mf_model->boxArray(),
                                                                         mf_model->DistributionMap(),
                                                                         mf_model->nComp(),
                                                                         amrex::IntVect::TheZeroVector() );
            }
        }

    }

    // Define the scalar data container
    if (m_scalar_type != FieldType::None) {

        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            const amrex::MultiFab* mf_model = a_WarpX->getFieldPointer(m_scalar_type,lev,0);
            m_scalar_vec[lev] = std::make_unique<amrex::MultiFab>( mf_model->boxArray(),
                                                                   mf_model->DistributionMap(),
                                                                   mf_model->nComp(),
                                                                   amrex::IntVect::TheZeroVector() );
        }

    }

    m_is_defined = true;

    // Define static member pointer to WarpX
    if (!m_warpx_ptr_defined) {
        m_WarpX = a_WarpX;
        m_warpx_ptr_defined = true;
    }
}

void WarpXSolverVec::Copy ( FieldType  a_array_type,
                            FieldType  a_scalar_type )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(),
        "WarpXSolverVec::Copy() called on undefined WarpXSolverVec");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_array_type==m_array_type &&
        a_scalar_type==m_scalar_type,
        "WarpXSolverVec::Copy() called with vecs of different types");

    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        if (m_array_type != FieldType::None) {
            for (int n = 0; n < 3; ++n) {
                const amrex::MultiFab* this_field = m_WarpX->getFieldPointer(m_array_type,lev,n);
                amrex::MultiFab::Copy( *m_array_vec[lev][n], *this_field, 0, 0, m_ncomp,
                                       amrex::IntVect::TheZeroVector() );
            }
        }
        if (m_scalar_type != FieldType::None) {
            const amrex::MultiFab* this_scalar = m_WarpX->getFieldPointer(m_scalar_type,lev,0);
            amrex::MultiFab::Copy( *m_scalar_vec[lev], *this_scalar, 0, 0, m_ncomp,
                                   amrex::IntVect::TheZeroVector() );
        }
    }
}

[[nodiscard]] amrex::Real WarpXSolverVec::dotProduct ( const WarpXSolverVec&  a_X ) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_X.getArrayVecType()==m_array_type &&
        a_X.getScalarVecType()==m_scalar_type,
        "WarpXSolverVec::dotProduct(X) called with solver vecs of different types");

    amrex::Real result = 0.0;
    const bool local = true;
    for (int lev = 0; lev < m_num_amr_levels; ++lev) {
        if (m_array_type != FieldType::None) {
            for (int n = 0; n < 3; ++n) {
                const amrex::iMultiFab* dotMask = m_WarpX->getFieldDotMaskPointer(m_array_type,lev,n);
                auto rtmp = amrex::MultiFab::Dot( *dotMask,
                                                  *m_array_vec[lev][n], 0,
                                                  *a_X.getArrayVec()[lev][n], 0, 1, 0, local);
                result += rtmp;
            }
        }
        if (m_scalar_type != FieldType::None) {
            const amrex::iMultiFab* dotMask = m_WarpX->getFieldDotMaskPointer(m_scalar_type,lev,0);
            auto rtmp = amrex::MultiFab::Dot( *dotMask,
                                              *m_scalar_vec[lev], 0,
                                              *a_X.getScalarVec()[lev], 0, 1, 0, local);
            result += rtmp;
        }
    }
    amrex::ParallelAllReduce::Sum(result, amrex::ParallelContext::CommunicatorSub());
    return result;
}
