/* This file is part of WarpX.
 *
 * Authors: Roelof Groenewald
 * License: BSD-3-Clause-LBNL
 */

#include "JFunctor.H"

#include "Particles/MultiParticleContainer.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

JFunctor::JFunctor (const int dir, int lev,
                   amrex::IntVect crse_ratio,
                   bool convertRZmodes2cartesian,
                   bool deposit_current, int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_dir(dir), m_lev(lev),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian),
      m_deposit_current(deposit_current)
{ }

void
JFunctor::operator() (amrex::MultiFab& mf_dst, int dcomp, const int /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();
    /** pointer to source multifab (can be multi-component) */
    amrex::MultiFab* m_mf_src = warpx.get_pointer_current_fp(m_lev, m_dir);

    // Deposit current if no solver or the electrostatic solver is being used
    if (m_deposit_current)
    {
        // allocate temporary multifab to deposit current density into
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp_temp;
        current_fp_temp.resize(1);

        const auto& current_fp_x = warpx.getcurrent_fp(m_lev,0);
        current_fp_temp[0][0] = std::make_unique<amrex::MultiFab>(
            current_fp_x, amrex::make_alias, 0, current_fp_x.nComp()
        );

        const auto& current_fp_y = warpx.getcurrent_fp(m_lev,1);
        current_fp_temp[0][1] = std::make_unique<amrex::MultiFab>(
            current_fp_y, amrex::make_alias, 0, current_fp_y.nComp()
        );
        const auto& current_fp_z = warpx.getcurrent_fp(m_lev,2);
        current_fp_temp[0][2] = std::make_unique<amrex::MultiFab>(
            current_fp_z, amrex::make_alias, 0, current_fp_z.nComp()
        );

        auto& mypc = warpx.GetPartContainer();
        mypc.DepositCurrent(current_fp_temp, warpx.getdt(m_lev), 0.0);
    }

    InterpolateMFForDiag(mf_dst, *m_mf_src, dcomp, warpx.DistributionMap(m_lev),
                         m_convertRZmodes2cartesian);
}
