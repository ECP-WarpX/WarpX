/* Copyright 2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "MaxField.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>

#include <iostream>
#include <cmath>


using namespace amrex;

// constructor
MaxField::MaxField (std::string rd_name)
: ReducedDiags{rd_name}
{

    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "MaxField reduced diagnostics does not work for RZ coordinate.");
#endif

    // read number of levels
    int nLevel = 0;
    ParmParse pp("amr");
    pp.query("max_level", nLevel);
    nLevel += 1;

    // resize data array
    m_data.resize(8*nLevel,0.0);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app);
            // write header row
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(3+8*lev) + "]";
                ofs << "max_Ex_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(4+8*lev) + "]";
                ofs << "max_Ey_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(5+8*lev) + "]";
                ofs << "max_Ez_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(6+8*lev) + "]";
                ofs << "max_|E|_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(7+8*lev) + "]";
                ofs << "max_Bx_lev"+std::to_string(lev)+" (T)";
                ofs << m_sep;
                ofs << "[" + std::to_string(8+8*lev) + "]";
                ofs << "max_By_lev"+std::to_string(lev)+" (T)";
                ofs << m_sep;
                ofs << "[" + std::to_string(9+8*lev) + "]";
                ofs << "max_Bz_lev"+std::to_string(lev)+" (T)";
                ofs << m_sep;
                ofs << "[" + std::to_string(10+8*lev) + "]";
                ofs << "max_|B|_lev"+std::to_string(lev)+" (T)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes maximum field values
void MaxField::ComputeDiags (int step)
{
    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get number of level
    auto nLevel = warpx.finestLevel() + 1;

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

        // get Maximums of E field components
        m_data[lev*8] = Ex.norm0();
        m_data[lev*8+1] = Ey.norm0();
        m_data[lev*8+2] = Ez.norm0();

        // get Maximums of B field components
        m_data[lev*8+4] = Bx.norm0();
        m_data[lev*8+5] = By.norm0();
        m_data[lev*8+6] = Bz.norm0();

        // Create temporary MultiFAB to be filled with |E| and |B| squared
        const int ncomp = 1;
        const int ngrow = 0;  // no ghost cells for temporary MultiFAB
        MultiFab mftemp(amrex::convert(Ex.boxArray(), IntVect{AMREX_D_DECL(0,0,0)}),
                                            Ex.DistributionMap(), ncomp, ngrow);
        // Temporary MultiFab is cell-centered so that it can be filled for any staggering
        // (possibly with the sum of components which do not have the same staggering).

        MFItInfo info;
        if (TilingIfNotGPU()) {
            info.EnableTiling();
        }
#ifdef _OPENMP
        info.SetDynamic(WarpX::do_dynamic_scheduling);
#pragma omp parallel
#endif
        // MFIter loop to fill temporary MultiFAB with |E| squared.
        for (MFIter mfi(mftemp, info); mfi.isValid(); ++mfi )
        {
            const Box& box = mfi.tilebox();

            const auto& arrEx = Ex[mfi].array();
            const auto& arrEy = Ey[mfi].array();
            const auto& arrEz = Ez[mfi].array();
            auto arrtemp = mftemp[mfi].array();

            amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i, int j, int k){
            arrtemp(i,j,k) = arrEx(i,j,k)*arrEx(i,j,k) + arrEy(i,j,k)*arrEy(i,j,k)
                               + arrEz(i,j,k)*arrEz(i,j,k);
        });
        }
        m_data[lev*8+3] = std::sqrt(mftemp.max(0));

        // MFIter loop to fill temporary MultiFAB with |B| squared.
        for (MFIter mfi(mftemp, info); mfi.isValid(); ++mfi )
        {
            const Box& box = mfi.tilebox();

            const auto& arrBx = Bx[mfi].array();
            const auto& arrBy = By[mfi].array();
            const auto& arrBz = Bz[mfi].array();
            auto arrtemp = mftemp[mfi].array();

            amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i, int j, int k){
            arrtemp(i,j,k) = arrBx(i,j,k)*arrBx(i,j,k) + arrBy(i,j,k)*arrBy(i,j,k)
                               + arrBz(i,j,k)*arrBz(i,j,k);
        });
        }
        m_data[lev*8+7] = std::sqrt(mftemp.max(0));
    }
    // end loop over refinement levels

    /* m_data now contains up-to-date values for:
     *  [max(Ex),max(Ey),max(Ez),max(|E|),
     *   max(Bx),max(By),max(Bz),max(|B|)] */

}
// end void MaxField::ComputeDiags
