/* This file is part of Warpx.
 *
 * Authors: Avigdor Veksler
 * License: BSD-3-Clause-LBNL
*/
#include "JeFunctor.H"

#include "WarpX.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

JeFunctor::JeFunctor (int dir, int lev, amrex::IntVect crse_ratio, int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_dir(dir), m_lev(lev)
{ }

void
JeFunctor::operator() (amrex::MultiFab& mf_dst, int dcomp, const int /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();
    auto& hybrid_pic_model = warpx.GetHybridPICModel();

    /** pointer to source1 (Ji) multifab */
    amrex::MultiFab* m_mf_j = warpx.get_pointer_current_fp(m_lev, m_dir);
    /** pointer to source2 (Jamp) multifab */
    amrex::MultiFab* m_mf_j_ampere = hybrid_pic_model.get_pointer_current_fp_ampere(m_lev, m_dir);
    /** pointer to source3 (Jext) multifab */
    amrex::MultiFab* m_mf_j_external = hybrid_pic_model.get_pointer_current_fp_external(m_lev, m_dir);

    // A Je multifab is generated to hold electron current. 
    amrex::MultiFab Je( m_mf_j->boxArray(), m_mf_j->DistributionMap(), 1, m_mf_j->nGrowVect() );

    // Je = Jamp - Ji - Jext
    amrex::MultiFab::LinComb(
        Je, 1, *m_mf_j_ampere, 0,
        -1, *m_mf_j, 0, 0, 1, Je.nGrowVect()
    );
    amrex::MultiFab::Subtract(Je, *m_mf_j_external, 0, 0, 1, Je.nGrowVect());

    InterpolateMFForDiag(mf_dst, Je, dcomp, warpx.DistributionMap(m_lev), false);
}
