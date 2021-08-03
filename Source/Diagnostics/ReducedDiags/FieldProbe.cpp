/* Copyright 2021 Lorenzo Giacomel
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldProbe.H"

#include "Utils/CoarsenIO.H"
#include "Utils/IntervalsParser.H"
#include "WarpX.H"

#include <AMReX_Array.H>
#include <AMReX_Config.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_Geometry.H>

#include <cmath>
#include <ostream>

using namespace amrex;

// constructor
FieldProbe::FieldProbe (std::string rd_name)
: ReducedDiags{rd_name}
{
    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldProbe reduced diagnostics does not work for RZ coordinate.");
#endif

    // read number of levels
    int nLevel = 0;
    ParmParse pp_amr("amr");
    pp_amr.query("max_level", nLevel);
    nLevel += 1;

    ParmParse pp_rd_name(rd_name);
    pp_rd_name.query("x_probe", x_probe);
    pp_rd_name.query("y_probe", y_probe);
    pp_rd_name.query("z_probe", z_probe);

    const bool raw_specified = pp_rd_name.query("raw_fields", raw_fields);
    if(!raw_specified) raw_fields = false;

    constexpr int noutputs = 8;  // probe Ex,Ey,Ez,|E|,Bx,By,Bz and |B|
    // resize data array
    m_data.resize(noutputs*nLevel, 0.0_rt);

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
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ex_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ey_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ez_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_|E|_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Bx_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_By_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Bz_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_|B|_lev" + std::to_string(lev) + " (T)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }
}
// end constructor

// function that computes maximum field values
void FieldProbe::ComputeDiags (int step)
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

        const amrex::Geometry& gm = warpx.Geom(lev);
        const auto prob_lo = gm.ProbLo();
        const auto cell_size = gm.CellSizeArray();

        const int i_probe = amrex::Math::floor((x_probe - prob_lo[0])/cell_size[0]);
        const int j_probe = amrex::Math::floor((y_probe - prob_lo[1])/cell_size[1]);
        const int k_probe = amrex::Math::floor((z_probe - prob_lo[2])/cell_size[2]);

        // get MultiFab data at lev
        const MultiFab & Ex = warpx.getEfield(lev,0);
        const MultiFab & Ey = warpx.getEfield(lev,1);
        const MultiFab & Ez = warpx.getEfield(lev,2);
        const MultiFab & Bx = warpx.getBfield(lev,0);
        const MultiFab & By = warpx.getBfield(lev,1);
        const MultiFab & Bz = warpx.getBfield(lev,2);

        constexpr int noutputs = 8; // probing Ex,Ey,Ez,|E|,Bx,By,Bz and |B|
        constexpr int index_Ex = 0;
        constexpr int index_Ey = 1;
        constexpr int index_Ez = 2;
        constexpr int index_absE = 3;
        constexpr int index_Bx = 4;
        constexpr int index_By = 5;
        constexpr int index_Bz = 6;
        constexpr int index_absB = 7;

        // General preparation of interpolation and reduction operations
        const GpuArray<int,3> cellCenteredtype{0,0,0};
        const GpuArray<int,3> reduction_coarsening_ratio{1,1,1};
        constexpr int reduction_comp = 0;

        // Prepare interpolation of field components to cell center
        // The arrays below store the index type (staggering) of each MultiFab, with the third
        // component set to zero in the two-dimensional case.
        auto Extype = amrex::GpuArray<int,3>{0,0,0};
        auto Eytype = amrex::GpuArray<int,3>{0,0,0};
        auto Eztype = amrex::GpuArray<int,3>{0,0,0};
        auto Bxtype = amrex::GpuArray<int,3>{0,0,0};
        auto Bytype = amrex::GpuArray<int,3>{0,0,0};
        auto Bztype = amrex::GpuArray<int,3>{0,0,0};
        for (int i = 0; i < AMREX_SPACEDIM; ++i){
            Extype[i] = Ex.ixType()[i];
            Eytype[i] = Ey.ixType()[i];
            Eztype[i] = Ez.ixType()[i];
            Bxtype[i] = Bx.ixType()[i];
            Bytype[i] = By.ixType()[i];
            Bztype[i] = Bz.ixType()[i];
        }

        amrex::Real hv_Ex, hv_Ey, hv_Ez, hv_E, hv_Bx, hv_By, hv_Bz, hv_B;
        int probe_proc = -1; //amrex::ParallelDescriptor::IOProcessorNumber();
        // MFIter loop to interpolate fields to cell center and get maximum values
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            const Box& box = mfi.validbox();

            if(box.contains(i_probe, j_probe, k_probe)) {
                // Make the box cell centered in preparation for the interpolation (and to avoid
                // including ghost cells in the calculation)
                const auto &arrEx = Ex[mfi].array();
                const auto &arrEy = Ey[mfi].array();
                const auto &arrEz = Ez[mfi].array();
                const auto &arrBx = Bx[mfi].array();
                const auto &arrBy = By[mfi].array();
                const auto &arrBz = Bz[mfi].array();

                //Retrieveing cell-centered fields
                const amrex::Real Ex_interp = CoarsenIO::Interp(arrEx, Extype, cellCenteredtype,
                                                                reduction_coarsening_ratio,
                                                                i_probe, j_probe, k_probe, reduction_comp);
                const amrex::Real Ey_interp = CoarsenIO::Interp(arrEy, Extype, cellCenteredtype,
                                                                reduction_coarsening_ratio,
                                                                i_probe, j_probe, k_probe, reduction_comp);
                const amrex::Real Ez_interp = CoarsenIO::Interp(arrEz, Extype, cellCenteredtype,
                                                                reduction_coarsening_ratio,
                                                                i_probe, j_probe, k_probe, reduction_comp);
                const amrex::Real Bx_interp = CoarsenIO::Interp(arrBx, Extype, cellCenteredtype,
                                                                reduction_coarsening_ratio,
                                                                i_probe, j_probe, k_probe, reduction_comp);
                const amrex::Real By_interp = CoarsenIO::Interp(arrBx, Extype, cellCenteredtype,
                                                                reduction_coarsening_ratio,
                                                                i_probe, j_probe, k_probe, reduction_comp);
                const amrex::Real Bz_interp = CoarsenIO::Interp(arrBz, Extype, cellCenteredtype,
                                                                reduction_coarsening_ratio,
                                                                i_probe, j_probe, k_probe, reduction_comp);

                // Either save the interpolated fields or the raw fields depending on the raw_fields flag
                hv_Ex = raw_fields ? arrEx(i_probe, j_probe, k_probe) : Ex_interp;
                hv_Ey = raw_fields ? arrEy(i_probe, j_probe, k_probe) : Ey_interp;
                hv_Ez = raw_fields ? arrEz(i_probe, j_probe, k_probe) : Ez_interp;
                hv_E = Ex_interp * Ex_interp + Ey_interp * Ey_interp + Ez_interp * Ez_interp;
                hv_Bx = raw_fields ? arrBx(i_probe, j_probe, k_probe) : Bx_interp;
                hv_By = raw_fields ? arrBy(i_probe, j_probe, k_probe) : By_interp;
                hv_Bz = raw_fields ? arrBz(i_probe, j_probe, k_probe) : Bz_interp;
                hv_B = Bx_interp * Bx_interp + By_interp * By_interp + Bz_interp * Bz_interp;

                probe_proc = ParallelDescriptor::MyProc();
            }
        }
        //All the processors have probe_proc = -1 except the one that contains the point, which
        //has probe_proc equal to a number >=0. In this way we communicate to all the processors
        //the rank of the processor which contains the point

        amrex::ParallelDescriptor::ReduceIntMax(probe_proc);
        Print()<<probe_proc<<std::endl;

        if(amrex::ParallelDescriptor::MyProc()==probe_proc){
            amrex::ParallelDescriptor::Send(&hv_Ex, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 0);
            amrex::ParallelDescriptor::Send(&hv_Ey, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 1);
            amrex::ParallelDescriptor::Send(&hv_Ez, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 2);
            amrex::ParallelDescriptor::Send(&hv_E, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 3);
            amrex::ParallelDescriptor::Send(&hv_Bx, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 4);
            amrex::ParallelDescriptor::Send(&hv_By, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 5);
            amrex::ParallelDescriptor::Send(&hv_Bz, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 6);
            amrex::ParallelDescriptor::Send(&hv_B, 1,
                                            amrex::ParallelDescriptor::IOProcessorNumber(), 7);
        }
        if(amrex::ParallelDescriptor::MyProc()==amrex::ParallelDescriptor::IOProcessorNumber()){
            amrex::ParallelDescriptor::Recv(&hv_Ex, 1, probe_proc, 0);
            amrex::ParallelDescriptor::Recv(&hv_Ey, 1, probe_proc, 1);
            amrex::ParallelDescriptor::Recv(&hv_Ez, 1, probe_proc, 2);
            amrex::ParallelDescriptor::Recv(&hv_E, 1, probe_proc, 3);
            amrex::ParallelDescriptor::Recv(&hv_Bx, 1, probe_proc, 4);
            amrex::ParallelDescriptor::Recv(&hv_By, 1, probe_proc, 5);
            amrex::ParallelDescriptor::Recv(&hv_Bz, 1, probe_proc, 6);
            amrex::ParallelDescriptor::Recv(&hv_B, 1, probe_proc, 7);

        }



        // Fill output array
        m_data[lev * noutputs + index_Ex] = hv_Ex;
        m_data[lev * noutputs + index_Ey] = hv_Ey;
        m_data[lev * noutputs + index_Ez] = hv_Ez;
        m_data[lev * noutputs + index_absE] = std::sqrt(hv_E);
        m_data[lev * noutputs + index_Bx] = hv_Bx;
        m_data[lev * noutputs + index_By] = hv_By;
        m_data[lev * noutputs + index_Bz] = hv_Bz;
        m_data[lev * noutputs + index_absB] = std::sqrt(hv_B);
    }
    // end loop over refinement levels

    /* m_data now contains up-to-date values for:
     *  [probe(Ex),probe(Ey),probe(Ez),probe(|E|),
     *   probe(Bx),probe(By),probe(Bz),probe(|B|)] */
}
// end void FieldProbe::ComputeDiags
