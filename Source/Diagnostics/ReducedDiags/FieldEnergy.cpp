/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldEnergy.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>

#include <iostream>
#include <cmath>


using namespace amrex;

// constructor
FieldEnergy::FieldEnergy (std::string rd_name)
: ReducedDiags{rd_name}
{

    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldEnergy reduced diagnostics does not work for RZ coordinate.");
#endif

    // read number of levels
    int nLevel = 0;
    ParmParse pp_amr("amr");
    pp_amr.query("max_level", nLevel);
    nLevel += 1;

    constexpr int noutputs = 3; // total energy, E-field energy and B-field energy
    // resize data array
    m_data.resize(noutputs*nLevel, 0.0_rt);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension, std::ofstream::out};
            // write header row
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            constexpr int shift_total = 3;
            constexpr int shift_E = 4;
            constexpr int shift_B = 5;
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_total+noutputs*lev) + "]";
                ofs << "total_lev"+std::to_string(lev)+"(J)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_E+noutputs*lev) + "]";
                ofs << "E_lev"+std::to_string(lev)+"(J)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_B+noutputs*lev) + "]";
                ofs << "B_lev"+std::to_string(lev)+"(J)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes field energy
void FieldEnergy::ComputeDiags (int step)
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
#if (AMREX_SPACEDIM == 2)
        auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif (AMREX_SPACEDIM == 3)
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
    // end loop over refinement levels

    /* m_data now contains up-to-date values for:
     *  [total field energy at level 0,
     *   electric field energy at level 0,
     *   magnetic field energy at level 0,
     *   total field energy at level 1,
     *   electric field energy at level 1,
     *   magnetic field energy at level 1,
     *   ......] */

}
// end void FieldEnergy::ComputeDiags
