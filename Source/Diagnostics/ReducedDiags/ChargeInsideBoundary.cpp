/* Copyright 2023 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ChargeInsideBoundary.H"

#include "Diagnostics/ReducedDiags/ReducedDiags.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <fstream>
#include <vector>

using namespace amrex;

// constructor
ChargeInsideBoundary::ChargeInsideBoundary (std::string rd_name)
: ReducedDiags{rd_name}
{
    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "ChargeInsideBoundary reduced diagnostics does not work for RZ coordinate.");
#endif

    // resize data array
    m_data.resize(1, 0.0_rt);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension, std::ofstream::out};
            // write header row
            int c = 0;
            ofs << "#";
            ofs << "[" << c++ << "]step()";
            ofs << m_sep;
            ofs << "[" << c++ << "]time(s)";
            ofs << m_sep;
            ofs << "[" << c++ << "]Charge (C)";
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }
}
// end constructor

// function that computes the charge inside a boundary
void ChargeInsideBoundary::ComputeDiags (int step)
{
    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // get number of level
    const auto nLevel = warpx.finestLevel() + 1;

    // loop over refinement levels
    for (int lev = 0; lev < nLevel; ++lev)
    {
        // get MultiFab data at lev
        const MultiFab & Ex = warpx.getEfield(lev,0);
        const MultiFab & Ey = warpx.getEfield(lev,1);
        const MultiFab & Ez = warpx.getEfield(lev,2);
        const MultiFab & Bx = warpx.getBfield(lev,0);
        const MultiFab & By = warpx.getBfield(lev,1);
        const MultiFab & Bz = warpx.getBfield(lev,2);

        // get cell size
        Geometry const & geom = warpx.Geom(lev);
#if defined(WARPX_DIM_1D_Z)
        auto dV = geom.CellSize(0);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined(WARPX_DIM_3D)
        auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

        // compute E squared
        Real const tmpEx = Ex.norm2(0,geom.periodicity());
        Real const tmpEy = Ey.norm2(0,geom.periodicity());
        Real const tmpEz = Ez.norm2(0,geom.periodicity());
        Real const Es = tmpEx*tmpEx + tmpEy*tmpEy + tmpEz*tmpEz;

        // compute B squared
        Real const tmpBx = Bx.norm2(0,geom.periodicity());
        Real const tmpBy = By.norm2(0,geom.periodicity());
        Real const tmpBz = Bz.norm2(0,geom.periodicity());
        Real const Bs = tmpBx*tmpBx + tmpBy*tmpBy + tmpBz*tmpBz;

        constexpr int noutputs = 3; // total energy, E-field energy and B-field energy
        constexpr int index_total = 0;
        constexpr int index_E = 1;
        constexpr int index_B = 2;

        // save data
        m_data[lev*noutputs+index_E] = 0.5_rt * Es * PhysConst::ep0 * dV;
        m_data[lev*noutputs+index_B] = 0.5_rt * Bs / PhysConst::mu0 * dV;
        m_data[lev*noutputs+index_total] = m_data[lev*noutputs+index_E] +
                                           m_data[lev*noutputs+index_B];
    }

}
// end void ChargeInsideBoundary::ComputeDiags
