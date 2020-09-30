/* Copyright 2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldMaximum.H"
#include "WarpX.H"
#include "Utils/CoarsenIO.H"

using namespace amrex;

// constructor
FieldMaximum::FieldMaximum (std::string rd_name)
: ReducedDiags{rd_name}
{

    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldMaximum reduced diagnostics does not work for RZ coordinate.");
#endif

    // read number of levels
    int nLevel = 0;
    ParmParse pp("amr");
    pp.query("max_level", nLevel);
    nLevel += 1;

    constexpr int noutputs = 8; // total energy, E-field energy and B-field energy
    // resize data array
    m_data.resize(noutputs*nLevel,0.0); // max of Ex,Ey,Ez,|E|,Bx,By,Bz and |B|

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
            constexpr int shift_Ex = 3;
            constexpr int shift_Ey = 4;
            constexpr int shift_Ez = 5;
            constexpr int shift_absE = 6;
            constexpr int shift_Bx = 7;
            constexpr int shift_By = 8;
            constexpr int shift_Bz = 9;
            constexpr int shift_absB = 10;
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_Ex+noutputs*lev) + "]";
                ofs << "max_Ex_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_Ey+noutputs*lev) + "]";
                ofs << "max_Ey_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_Ez+noutputs*lev) + "]";
                ofs << "max_Ez_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_absE+noutputs*lev) + "]";
                ofs << "max_|E|_lev"+std::to_string(lev)+" (V/m)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_Bx+noutputs*lev) + "]";
                ofs << "max_Bx_lev"+std::to_string(lev)+" (T)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_By+noutputs*lev) + "]";
                ofs << "max_By_lev"+std::to_string(lev)+" (T)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_Bz+noutputs*lev) + "]";
                ofs << "max_Bz_lev"+std::to_string(lev)+" (T)";
                ofs << m_sep;
                ofs << "[" + std::to_string(shift_absB+noutputs*lev) + "]";
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
void FieldMaximum::ComputeDiags (int step)
{
    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get WarpX class object
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

        // get Maximums of E field components
        m_data[lev*noutputs+index_Ex] = Ex.norm0();
        m_data[lev*noutputs+index_Ey] = Ey.norm0();
        m_data[lev*noutputs+index_Ez] = Ez.norm0();

        // get Maximums of B field components
        m_data[lev*noutputs+index_Bx] = Bx.norm0();
        m_data[lev*noutputs+index_By] = By.norm0();
        m_data[lev*noutputs+index_Bz] = Bz.norm0();

        // General preparation of interpolation and reduction to compute the maximum value of
        // |E| and |B|
        const GpuArray<int,3> cellCenteredtype{0,0,0};
        const GpuArray<int,3> reduction_coarsening_ratio{1,1,1};
        constexpr int reduction_comp = 0;

        // Prepare reduction for maximum value of |E|
        ReduceOps<ReduceOpMax> reduceE_op;
        ReduceData<Real> reduceE_data(reduceE_op);
        using ReduceTuple = typename decltype(reduceE_data)::Type;

        // Prepare interpolation of E field components to cell center
        GpuArray<int,3> Extype;
        GpuArray<int,3> Eytype;
        GpuArray<int,3> Eztype;
        for (int i = 0; i < AMREX_SPACEDIM; ++i){
            Extype[i] = Ex.ixType()[i];
            Eytype[i] = Ey.ixType()[i];
            Eztype[i] = Ez.ixType()[i];
        }
#if   (AMREX_SPACEDIM == 2)
        Extype[2] = 0;
        Eytype[2] = 0;
        Eztype[2] = 0;
#endif

        // MFIter loop to interpolate E field components to cell center and get maximum value
        // of |E|
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Make the box cell centered to avoid including ghost cells in the calculation
            const Box& box = enclosedCells(mfi.nodaltilebox());
            const auto& arrEx = Ex[mfi].array();
            const auto& arrEy = Ey[mfi].array();
            const auto& arrEz = Ez[mfi].array();
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
        }

        Real hv_E = amrex::get<0>(reduceE_data.value()); // highest value of |E|**2
        // MPI reduce
        ParallelDescriptor::ReduceRealMax(hv_E);

        m_data[lev*noutputs+index_absE] = std::sqrt(hv_E);

        // Prepare reduction for maximum value of |B|
        ReduceOps<ReduceOpMax> reduceB_op;
        ReduceData<Real> reduceB_data(reduceB_op);
        using ReduceTuple = typename decltype(reduceB_data)::Type;

        // Prepare interpolation of B field components to cell center
        GpuArray<int,3> Bxtype;
        GpuArray<int,3> Bytype;
        GpuArray<int,3> Bztype;
        for (int i = 0; i < AMREX_SPACEDIM; ++i){
            Bxtype[i] = Bx.ixType()[i];
            Bytype[i] = By.ixType()[i];
            Bztype[i] = Bz.ixType()[i];
        }
#if   (AMREX_SPACEDIM == 2)
        Bxtype[2] = 0;
        Bytype[2] = 0;
        Bztype[2] = 0;
#endif

        // MFIter loop to interpolate B field components to cell center and get maximum value
        // of |B|
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Make the box cell centered to avoid including ghost cells in the calculation
            const Box& box = enclosedCells(mfi.nodaltilebox());
            const auto& arrBx = Bx[mfi].array();
            const auto& arrBy = By[mfi].array();
            const auto& arrBz = Bz[mfi].array();

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

        Real hv_B = amrex::get<0>(reduceB_data.value()); // highest value of |B|**2
        // MPI reduce
        ParallelDescriptor::ReduceRealMax(hv_B);

        m_data[lev*noutputs+index_absB] = std::sqrt(hv_B);
    }
    // end loop over refinement levels

    /* m_data now contains up-to-date values for:
     *  [max(Ex),max(Ey),max(Ez),max(|E|),
     *   max(Bx),max(By),max(Bz),max(|B|)] */

}
// end void FieldMaximum::ComputeDiags
