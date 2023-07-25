/* This file is part of WarpX.
 *
 * Authors: Roelof Groenewald
 * License: BSD-3-Clause-LBNL
 */

#include "JFunctor.H"

#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

JFunctor::JFunctor(const int dir, int lev,
                   amrex::IntVect crse_ratio,
                   bool convertRZmodes2cartesian, int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_dir(dir), m_lev(lev),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{ }

void
JFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp, const int /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();
    /** pointer to source multifab (can be multi-component) */
    amrex::MultiFab* m_mf_src = warpx.get_pointer_current_fp(m_lev, m_dir);

    // Deposit current if no solver or the electrostatic solver is being used
    if ( WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::None &&
        WarpX::electrostatic_solver_id != ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic
    ) {
        // allocate temporary multifab to deposit current density into
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp_temp;
        current_fp_temp.resize(1);

        auto& current_fp_x = warpx.getcurrent_fp(m_lev,0);
        current_fp_temp[0][0] = std::make_unique<amrex::MultiFab>(
            current_fp_x.boxArray(), current_fp_x.DistributionMap(),
            current_fp_x.nComp(), current_fp_x.nGrow()
        );
        auto& current_fp_y = warpx.getcurrent_fp(m_lev,1);
        current_fp_temp[0][1] = std::make_unique<amrex::MultiFab>(
            current_fp_y.boxArray(), current_fp_y.DistributionMap(),
            current_fp_y.nComp(), current_fp_y.nGrow()
        );
        auto& current_fp_z = warpx.getcurrent_fp(m_lev,2);
        current_fp_temp[0][2] = std::make_unique<amrex::MultiFab>(
            current_fp_z.boxArray(), current_fp_z.DistributionMap(),
            current_fp_z.nComp(), current_fp_z.nGrow()
        );

        auto& mypc = warpx.GetPartContainer();
        mypc.DepositCurrent(current_fp_temp, warpx.getdt(m_lev), 0.0);

        // copy the deposited current into the src mf
        amrex::MultiFab::Copy(*m_mf_src, *current_fp_temp[0][m_dir], 0, 0, 1, m_mf_src->nGrowVect());
    }

#ifdef WARPX_DIM_RZ
    if (m_convertRZmodes2cartesian) {
        // In cylindrical geometry, sum real part of all modes of m_mf_src in
        // temporary multifab mf_dst_stag, and cell-center it to mf_dst.
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into 1 single component");
        amrex::MultiFab mf_dst_stag(m_mf_src->boxArray(), warpx.DistributionMap(m_lev), 1, m_mf_src->nGrowVect());
        // Mode 0
        amrex::MultiFab::Copy(mf_dst_stag, *m_mf_src, 0, 0, 1, m_mf_src->nGrowVect());
        for (int ic=1 ; ic < m_mf_src->nComp() ; ic += 2) {
            // All modes > 0
            amrex::MultiFab::Add(mf_dst_stag, *m_mf_src, ic, 0, 1, m_mf_src->nGrowVect());
        }
        ablastr::coarsen::sample::Coarsen( mf_dst, mf_dst_stag, dcomp, 0, nComp(), 0,  m_crse_ratio);
    } else {
        ablastr::coarsen::sample::Coarsen( mf_dst, *m_mf_src, dcomp, 0, nComp(), 0, m_crse_ratio);
    }
#else
    // In cartesian geometry, coarsen and interpolate from simulation MultiFab, m_mf_src,
    // to output diagnostic MultiFab, mf_dst.
    ablastr::coarsen::sample::Coarsen(mf_dst, *m_mf_src, dcomp, 0, nComp(), mf_dst.nGrowVect(), m_crse_ratio);
    amrex::ignore_unused(m_lev, m_convertRZmodes2cartesian);
#endif
}
