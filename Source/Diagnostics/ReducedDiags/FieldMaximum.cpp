/* Copyright 2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldMaximum.H"

#include "Utils/CoarsenIO.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

using namespace amrex;

// constructor
FieldMaximum::FieldMaximum (std::string rd_name)
: ReducedDiags{rd_name}
{
    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldMaximum reduced diagnostics does not work for RZ coordinate.");
#endif

    // read number of levels
    int nLevel = 0;
    ParmParse pp_amr("amr");
    pp_amr.query("max_level", nLevel);
    nLevel += 1;

    constexpr int noutputs = 8;  // max of Ex,Ey,Ez,|E|,Bx,By,Bz and |B|
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
                ofs << "[" << c++ << "]max_Ex_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]max_Ey_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]max_Ez_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]max_|E|_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]max_Bx_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]max_By_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]max_Bz_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]max_|B|_lev" + std::to_string(lev) + " (T)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }
}
// end constructor

// function that computes maximum field values
void FieldMaximum::ComputeDiags (int step)
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

        constexpr int noutputs = 8; // max of Ex,Ey,Ez,|E|,Bx,By,Bz and |B|
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

        ReduceOps<ReduceOpMax> reduceEx_op;
        ReduceOps<ReduceOpMax> reduceEy_op;
        ReduceOps<ReduceOpMax> reduceEz_op;
        ReduceOps<ReduceOpMax> reduceBx_op;
        ReduceOps<ReduceOpMax> reduceBy_op;
        ReduceOps<ReduceOpMax> reduceBz_op;
        ReduceOps<ReduceOpMax> reduceE_op;
        ReduceOps<ReduceOpMax> reduceB_op;

        ReduceData<Real> reduceEx_data(reduceEx_op);
        ReduceData<Real> reduceEy_data(reduceEy_op);
        ReduceData<Real> reduceEz_data(reduceEz_op);
        ReduceData<Real> reduceBx_data(reduceBx_op);
        ReduceData<Real> reduceBy_data(reduceBy_op);
        ReduceData<Real> reduceBz_data(reduceBz_op);
        ReduceData<Real> reduceE_data(reduceE_op);
        ReduceData<Real> reduceB_data(reduceB_op);

        using ReduceTuple = typename decltype(reduceEx_data)::Type;

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

        // MFIter loop to interpolate fields to cell center and get maximum values
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Make the box cell centered in preparation for the interpolation (and to avoid
            // including ghost cells in the calculation)
            const Box& box = enclosedCells(mfi.nodaltilebox());
            const auto& arrEx = Ex[mfi].array();
            const auto& arrEy = Ey[mfi].array();
            const auto& arrEz = Ez[mfi].array();
            const auto& arrBx = Bx[mfi].array();
            const auto& arrBy = By[mfi].array();
            const auto& arrBz = Bz[mfi].array();

            reduceEx_op.eval(box, reduceEx_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real Ex_interp = CoarsenIO::Interp(arrEx, Extype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return amrex::Math::abs(Ex_interp);
            });
            reduceEy_op.eval(box, reduceEy_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real Ey_interp = CoarsenIO::Interp(arrEy, Eytype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return amrex::Math::abs(Ey_interp);
            });
            reduceEz_op.eval(box, reduceEz_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real Ez_interp = CoarsenIO::Interp(arrEz, Eztype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return amrex::Math::abs(Ez_interp);
            });
            reduceBx_op.eval(box, reduceBx_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real Bx_interp = CoarsenIO::Interp(arrBx, Bxtype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return amrex::Math::abs(Bx_interp);
            });
            reduceBy_op.eval(box, reduceBy_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real By_interp = CoarsenIO::Interp(arrBy, Bytype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return amrex::Math::abs(By_interp);
            });
            reduceBz_op.eval(box, reduceBz_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real Bz_interp = CoarsenIO::Interp(arrBz, Bztype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return amrex::Math::abs(Bz_interp);
            });
            reduceE_op.eval(box, reduceE_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real Ex_interp = CoarsenIO::Interp(arrEx, Extype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                const Real Ey_interp = CoarsenIO::Interp(arrEy, Eytype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                const Real Ez_interp = CoarsenIO::Interp(arrEz, Eztype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return Ex_interp*Ex_interp + Ey_interp*Ey_interp + Ez_interp*Ez_interp;
            });
            reduceB_op.eval(box, reduceB_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real Bx_interp = CoarsenIO::Interp(arrBx, Bxtype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                const Real By_interp = CoarsenIO::Interp(arrBy, Bytype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                const Real Bz_interp = CoarsenIO::Interp(arrBz, Bztype, cellCenteredtype,
                                        reduction_coarsening_ratio, i, j, k, reduction_comp);
                return Bx_interp*Bx_interp + By_interp*By_interp + Bz_interp*Bz_interp;
            });
        }

        Real hv_Ex = amrex::get<0>(reduceEx_data.value()); // highest value of |Ex|
        Real hv_Ey = amrex::get<0>(reduceEy_data.value()); // highest value of |Ey|
        Real hv_Ez = amrex::get<0>(reduceEz_data.value()); // highest value of |Ez|
        Real hv_Bx = amrex::get<0>(reduceBx_data.value()); // highest value of |Bx|
        Real hv_By = amrex::get<0>(reduceBy_data.value()); // highest value of |By|
        Real hv_Bz = amrex::get<0>(reduceBz_data.value()); // highest value of |Bz|
        Real hv_E = amrex::get<0>(reduceE_data.value()); // highest value of |E|**2
        Real hv_B = amrex::get<0>(reduceB_data.value()); // highest value of |B|**2

        // MPI reduce
        ParallelDescriptor::ReduceRealMax(hv_Ex);
        ParallelDescriptor::ReduceRealMax(hv_Ey);
        ParallelDescriptor::ReduceRealMax(hv_Ez);
        ParallelDescriptor::ReduceRealMax(hv_Bx);
        ParallelDescriptor::ReduceRealMax(hv_By);
        ParallelDescriptor::ReduceRealMax(hv_Bz);
        ParallelDescriptor::ReduceRealMax(hv_E);
        ParallelDescriptor::ReduceRealMax(hv_B);

        // Fill output array
        m_data[lev*noutputs+index_Ex] = hv_Ex;
        m_data[lev*noutputs+index_Ey] = hv_Ey;
        m_data[lev*noutputs+index_Ez] = hv_Ez;
        m_data[lev*noutputs+index_Bx] = hv_Bx;
        m_data[lev*noutputs+index_By] = hv_By;
        m_data[lev*noutputs+index_Bz] = hv_Bz;
        m_data[lev*noutputs+index_absE] = std::sqrt(hv_E);
        m_data[lev*noutputs+index_absB] = std::sqrt(hv_B);
    }
    // end loop over refinement levels

    /* m_data now contains up-to-date values for:
     *  [max(Ex),max(Ey),max(Ez),max(|E|),
     *   max(Bx),max(By),max(Bz),max(|B|)] */
}
// end void FieldMaximum::ComputeDiags
