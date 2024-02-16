/* This file is part of Warpx.
 *
 * Authors: Avigdor Veksler
 * License: BSD-3-Clause-LBNL
*/
#include "JdispFunctor.H"

#include "WarpX.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

JdispFunctor::JdispFunctor (int dir, int lev,
        amrex::IntVect crse_ratio, bool convertRZmodes2cartesian, int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_dir(dir), m_lev(lev),
    m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{ }

void
JdispFunctor::operator() (amrex::MultiFab& mf_dst, int dcomp, const int /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();
    auto* hybrid_pic_model = warpx.get_pointer_HybridPICModel();

    /** pointer to source1 (total simulation current J) multifab */
    amrex::MultiFab* m_mf_j = warpx.get_pointer_current_fp(m_lev, m_dir);
    amrex::MultiFab* m_mf_curlB;
    amrex::MultiFab* m_mf_j_external;
    if (hybrid_pic_model) {
        /** pointer to source2 (current calculated from Ampere's Law, Jamp) multifab */
        m_mf_curlB = hybrid_pic_model->get_pointer_current_fp_ampere(m_lev, m_dir);
        /** pointer to source3 (external currents, Jext) multifab */
        m_mf_j_external = hybrid_pic_model->get_pointer_current_fp_external(m_lev, m_dir);
    } else {
        // To finish this implementation, we need to implement a method to
        // calculate (∇ x B).
        WARPX_ABORT_WITH_MESSAGE(
            "Displacement current diagnostic is only implemented for the HybridPICModel.");
            
        // Skeleton for future implementation for solvers other than HybridPIC.
        // Get curlB multifab

        // Multiply J multifab by mu0
        m_mf_curlB->mult(1.0/PhysConsts::mu0) */
    }
    // A Jdisp multifab is generated to hold displacement current.
    amrex::MultiFab Jdisp( m_mf_j->boxArray(), m_mf_j->DistributionMap(), 1, m_mf_j->nGrowVect() );

    // J_displacement = curl x B / mu0 - J
    amrex::MultiFab::LinComb(
        Jdisp, 1, *m_mf_curlB, 0,
        -1, *m_mf_j, 0, 0, 1, Jdisp.nGrowVect()
    );
    
    // subtract external current in the hybrid-PIC case
    if (hybrid_pic_model) {
        amrex::MultiFab::Subtract(Jdisp, *m_mf_j_external, 0, 0, 1, Jdisp.nGrowVect());
    } 

    InterpolateMFForDiag(mf_dst, Jdisp, dcomp, warpx.DistributionMap(m_lev),
                         m_convertRZmodes2cartesian);
}
