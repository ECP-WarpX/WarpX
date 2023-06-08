/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDREDUCTION_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDREDUCTION_H_

#include "ReducedDiags.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>

#include <memory>
#include <string>
#include <type_traits>
#include <vector>


/**
 * This class contains a function that computes an arbitrary reduction of the fields. The function
 * used in the reduction is defined by an input file parser expression and the reduction operation
 * can be either Maximum, Minimum, or Integral (Sum multiplied by cell volume).
 */
class FieldReduction : public ReducedDiags
{
public:

    /**
     * constructor
     * @param[in] rd_name reduced diags names
     */
    FieldReduction(std::string rd_name);

    /**
     * This function is called at every time step, and if necessary calls the templated function
     * ComputeFieldReduction(), which does the actual reduction computation.
     *
     * @param[in] step the timestep
     */
    virtual void ComputeDiags(int step) override final;

private:
    /// Parser to read expression to be reduced from the input file.
    /// 9 elements are x, y, z, Ex, Ey, Ez, Bx, By, Bz
    static constexpr int m_nvars = 9;
    std::unique_ptr<amrex::Parser> m_parser;

    // Type of reduction (e.g. Maximum, Minimum or Sum)
    int m_reduction_type;

public:

    /**
     * This function does the actual reduction computation. The fields are first interpolated on
     * the cell centers and the reduction operation is then performed using amrex::ReduceOps.
     *
     * \tparam ReduceOp the type of reduction that is performed. This is typically
     *         amrex::ReduceOpMax, amrex::ReduceOpMin or amrex::ReduceOpSum.
     */
    template<typename ReduceOp>
    void ComputeFieldReduction()
    {
        using namespace amrex::literals;

        // get a reference to WarpX instance
        auto & warpx = WarpX::GetInstance();

        constexpr int lev = 0; // This reduced diag currently does not work with mesh refinement

        amrex::Geometry const & geom = warpx.Geom(lev);
        const amrex::RealBox& real_box = geom.ProbDomain();
        const auto dx = geom.CellSizeArray();

        // get MultiFab data
        const amrex::MultiFab & Ex = warpx.getEfield(lev,0);
        const amrex::MultiFab & Ey = warpx.getEfield(lev,1);
        const amrex::MultiFab & Ez = warpx.getEfield(lev,2);
        const amrex::MultiFab & Bx = warpx.getBfield(lev,0);
        const amrex::MultiFab & By = warpx.getBfield(lev,1);
        const amrex::MultiFab & Bz = warpx.getBfield(lev,2);

        // General preparation of interpolation and reduction operations
        const amrex::GpuArray<int,3> cellCenteredtype{0,0,0};
        const amrex::GpuArray<int,3> reduction_coarsening_ratio{1,1,1};
        constexpr int reduction_comp = 0;

        amrex::ReduceOps<ReduceOp> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

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

        // get parser
        auto reduction_function_parser = m_parser->compile<m_nvars>();

        // MFIter loop to interpolate fields to cell center and perform reduction
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( amrex::MFIter mfi(Ex, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Make the box cell centered in preparation for the interpolation (and to avoid
            // including ghost cells in the calculation)
            const amrex::Box& box = enclosedCells(mfi.nodaltilebox());
            const auto& arrEx = Ex[mfi].array();
            const auto& arrEy = Ey[mfi].array();
            const auto& arrEz = Ez[mfi].array();
            const auto& arrBx = Bx[mfi].array();
            const auto& arrBy = By[mfi].array();
            const auto& arrBz = Bz[mfi].array();

            reduce_op.eval(box, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                 // 0.5 is here because position are computed on the cell centers

#if defined(WARPX_DIM_1D_Z)
                const amrex::Real x = 0._rt;
                const amrex::Real y = 0._rt;
                const amrex::Real z = (k + 0.5_rt)*dx[0] + real_box.lo(0);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const amrex::Real x = (i + 0.5_rt)*dx[0] + real_box.lo(0);
                const amrex::Real y = 0._rt;
                const amrex::Real z = (j + 0.5_rt)*dx[1] + real_box.lo(1);
#else
                const amrex::Real x = (i + 0.5_rt)*dx[0] + real_box.lo(0);
                const amrex::Real y = (j + 0.5_rt)*dx[1] + real_box.lo(1);
                const amrex::Real z = (k + 0.5_rt)*dx[2] + real_box.lo(2);
#endif
                const amrex::Real Ex_interp = ablastr::coarsen::sample::Interp(arrEx, Extype, cellCenteredtype,
                                                                               reduction_coarsening_ratio, i, j, k, reduction_comp);
                const amrex::Real Ey_interp = ablastr::coarsen::sample::Interp(arrEy, Eytype, cellCenteredtype,
                                                                               reduction_coarsening_ratio, i, j, k, reduction_comp);
                const amrex::Real Ez_interp = ablastr::coarsen::sample::Interp(arrEz, Eztype, cellCenteredtype,
                                                                               reduction_coarsening_ratio, i, j, k, reduction_comp);
                const amrex::Real Bx_interp = ablastr::coarsen::sample::Interp(arrBx, Bxtype, cellCenteredtype,
                                                                               reduction_coarsening_ratio, i, j, k, reduction_comp);
                const amrex::Real By_interp = ablastr::coarsen::sample::Interp(arrBy, Bytype, cellCenteredtype,
                                                                               reduction_coarsening_ratio, i, j, k, reduction_comp);
                const amrex::Real Bz_interp = ablastr::coarsen::sample::Interp(arrBz, Bztype, cellCenteredtype,
                                                                               reduction_coarsening_ratio, i, j, k, reduction_comp);
                return reduction_function_parser(x, y, z, Ex_interp, Ey_interp, Ez_interp,
                                                          Bx_interp, By_interp, Bz_interp);
            });
        }

        amrex::Real reduce_value = amrex::get<0>(reduce_data.value());

        // MPI reduce
        if (std::is_same<ReduceOp, amrex::ReduceOpMax>::value)
        {
            amrex::ParallelDescriptor::ReduceRealMax(reduce_value);
        }
        if (std::is_same<ReduceOp, amrex::ReduceOpMin>::value)
        {
            amrex::ParallelDescriptor::ReduceRealMin(reduce_value);
        }
        if (std::is_same<ReduceOp, amrex::ReduceOpSum>::value)
        {
            amrex::ParallelDescriptor::ReduceRealSum(reduce_value);
        // If reduction operation is a sum, multiply the value by the cell volume so that the
        // result is the integral of the function over the simulation domain.
#if defined(WARPX_DIM_1D_Z)
            reduce_value *= dx[0];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            reduce_value *= dx[0]*dx[1];
#else
            reduce_value *= dx[0]*dx[1]*dx[2];
#endif
        }

        // Fill output array
        m_data[0] = reduce_value;

        // m_data now contains an up-to-date value of the reduced field quantity
    }

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDREDUCTION_H_
