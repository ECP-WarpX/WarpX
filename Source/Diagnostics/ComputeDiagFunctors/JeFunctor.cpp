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

JeFunctor::JeFunctor (const int dir, int lev,
                     amrex::IntVect crse_ratio,
                     bool deposit_current, int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_dir(dir), m_lev(lev),
      m_deposit_current(deposit_current)
{ }

void
JeFunctor::operator() (amrex::MultiFab& mf_dst, int dcomp, const in t /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();
    /** pointer to source1 multifab */
    amrex::MultiFab* m_mf_j = warpx.get_pointer_current_fp(m_lev, m_dir);
    /** pointer to source2 multifab */
    amrex::MultiFab* m_mf_j_ampere = warpx.get_pointer_current_fp_ampere(m_lev, m_dir);

    // Deposit current if no solver is being used
    if (m_deposit_current)
    {
        // allocate temporary multifab to deposit current density into
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp_temp;
        current_fp_temp.resize(1);

        auto& current_fp_x = warpx.getcurrent_fp(m_lev,0);
        current_fp_temp[0][0] = std::make_unique<amrex::MultiFab>(
            current_fp_x, amrex::make_alias, 0, current_fp_x.nComp()
        );

        auto& current_fp_y = warpx.getcurrent_fp(m_lev,1);
        current_fp_temp[0][1] = std::make_unique<amrex::MultiFab>(
            current_fp_y, amrex::make_alias, 0, current_fp_y.nComp()
        );
        auto& current_fp_z = warpx.getcurrent_fp(m_lev,2);
        current_fp_temp[0][2] = std::make_unique<amrex::MultiFab>(
            current_fp_z, amrex::make_alias, 0, current_fp_z.nComp()
        );

        auto& mypc = warpx.GetPartContainer();
        mypc.DepositCurrent(current_fp_temp, warpx.getdt(m_lev), 0.0);
    }
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed.
    constexpr int ng = 1;
    // A cell-centered Je multifab spanning the entire domain is generated
    // to hold electron current. It has 1 component as the RZ Hybrid model
    // only supports the m=0 mode.
    amrex::MultiFab Je( warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 1, ng );

    // Subtract ampere current from ion current and store result Je multifab.
    MultiFab::LinComb(
        Je[m_lev], 1, (*m_mf_j)[m_lev], 0,
        -1, (*m_mf_j_ampere)[lev], 0, 0, 1, ng
    );

    InterpolateMFForDiag(mf_dst, Je, dcomp, warpx.DistributionMap(m_lev));
}