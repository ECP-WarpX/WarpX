/* Copyright 2019-2020
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldMomentum.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "Utils/CoarsenIO.H"

#include <AMReX_REAL.H>

#include <iostream>
#include <cmath>

using namespace amrex;

FieldMomentum::FieldMomentum (std::string rd_name)
    : ReducedDiags{rd_name}
{
    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
        amrex::Abort("FieldMomentum reduced diagnostics not implemented in RZ geometry");
#endif

    // read number of levels
    int nLevel = 0;
    amrex::ParmParse pp_amr("amr");
    pp_amr.query("max_level", nLevel);
    nLevel += 1;

    // Resize data array
    m_data.resize(nLevel*3, 0.0_rt);

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if (m_IsNotRestart)
        {
            // Open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension, std::ofstream::out};

            int c = 0;

            // Write header row
            ofs << "#";
            ofs << "[" << c++ << "]step()";
            ofs << m_sep;
            ofs << "[" << c++ << "]time(s)";

            // Loop over refinement levels
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" << c++ << "]";
                ofs << "momentum_x_lev" << lev << "(kg*m/s)";
                ofs << m_sep;
                ofs << "[" << c++ << "]";
                ofs << "momentum_y_lev" << lev << "(kg*m/s)";
                ofs << m_sep;
                ofs << "[" << c++ << "]";
                ofs << "momentum_z_lev" << lev << "(kg*m/s)";
            }

            ofs << std::endl;
            ofs.close();
        }
    }
}

void FieldMomentum::ComputeDiags (int step)
{
    // Check if the diags should be done
    if (m_intervals.contains(step+1) == false)
    {
        return;
    }

    // Get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // Get number of refinement levels
    const auto nLevel = warpx.finestLevel() + 1;

    // Loop over refinement levels
    for (int lev = 0; lev < nLevel; ++lev)
    {
        // Get MultiFab data at given refinement level
        const amrex::MultiFab & Ex = warpx.getEfield(lev, 0);
        const amrex::MultiFab & Ey = warpx.getEfield(lev, 1);
        const amrex::MultiFab & Ez = warpx.getEfield(lev, 2);
        const amrex::MultiFab & Bx = warpx.getBfield(lev, 0);
        const amrex::MultiFab & By = warpx.getBfield(lev, 1);
        const amrex::MultiFab & Bz = warpx.getBfield(lev, 2);

        // Allocate cell-centered MultiFabs
        amrex::MultiFab Ex_cc(amrex::convert(Ex.boxArray(), amrex::IntVect::TheZeroVector()),
                              warpx.DistributionMap(lev), 1, 0);

        amrex::MultiFab Ey_cc(amrex::convert(Ey.boxArray(), amrex::IntVect::TheZeroVector()),
                              warpx.DistributionMap(lev), 1, 0);

        amrex::MultiFab Ez_cc(amrex::convert(Ez.boxArray(), amrex::IntVect::TheZeroVector()),
                              warpx.DistributionMap(lev), 1, 0);

        amrex::MultiFab Bx_cc(amrex::convert(Bx.boxArray(), amrex::IntVect::TheZeroVector()),
                              warpx.DistributionMap(lev), 1, 0);

        amrex::MultiFab By_cc(amrex::convert(By.boxArray(), amrex::IntVect::TheZeroVector()),
                              warpx.DistributionMap(lev), 1, 0);

        amrex::MultiFab Bz_cc(amrex::convert(Bz.boxArray(), amrex::IntVect::TheZeroVector()),
                              warpx.DistributionMap(lev), 1, 0);

        // Interpolate to cell-centered MultiFabs
        CoarsenIO::Coarsen(Ex_cc, Ex, 0, 0, 1, 0);
        CoarsenIO::Coarsen(Ey_cc, Ey, 0, 0, 1, 0);
        CoarsenIO::Coarsen(Ez_cc, Ez, 0, 0, 1, 0);
        CoarsenIO::Coarsen(Bx_cc, Bx, 0, 0, 1, 0);
        CoarsenIO::Coarsen(By_cc, By, 0, 0, 1, 0);
        CoarsenIO::Coarsen(Bz_cc, Bz, 0, 0, 1, 0);

        // Get cell size
        amrex::Geometry const & geom = warpx.Geom(lev);
#if   (AMREX_SPACEDIM == 2)
        auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif (AMREX_SPACEDIM == 3)
        auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

        // Compute E x B (including sum over cells)
        const bool local = false;
        const amrex::Real ExB_x = amrex::MultiFab::Dot(Ey_cc, 0, Bz_cc, 0, 1, 0, local)
                                - amrex::MultiFab::Dot(Ez_cc, 0, By_cc, 0, 1, 0, local);
        const amrex::Real ExB_y = amrex::MultiFab::Dot(Ez_cc, 0, Bx_cc, 0, 1, 0, local)
                                - amrex::MultiFab::Dot(Ex_cc, 0, Bz_cc, 0, 1, 0, local);
        const amrex::Real ExB_z = amrex::MultiFab::Dot(Ex_cc, 0, By_cc, 0, 1, 0, local)
                                - amrex::MultiFab::Dot(Ey_cc, 0, Bx_cc, 0, 1, 0, local);

        // Save data (offset: 3 values for each refinement level)
        const int offset = lev*3;
        m_data[offset+0] = PhysConst::ep0 * ExB_x * dV;
        m_data[offset+1] = PhysConst::ep0 * ExB_y * dV;
        m_data[offset+2] = PhysConst::ep0 * ExB_z * dV;
    }
}
