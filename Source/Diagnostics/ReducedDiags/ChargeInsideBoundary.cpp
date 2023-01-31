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

#include <AMReX_GpuAtomic.H>
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
    // Only 3D is working for now
#if !(defined WARPX_DIM_3D)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "ChargeInsideBoundary reduced diagnostics only works in 3D");
#endif

#if !(defined AMREX_USE_EB)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "ChargeInsideBoundary reduced diagnostics only works when compiling with EB support");
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
#if ((defined WARPX_DIM_3D) && (defined AMREX_USE_EB))
    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // Only compute the integral on level 0
    int const lev = 0;

    // get MultiFab data at lev
    const MultiFab & Ex = warpx.getEfield(lev,0);
    const MultiFab & Ey = warpx.getEfield(lev,1);
    const MultiFab & Ez = warpx.getEfield(lev,2);

    // get EB structures
    amrex::EBFArrayBoxFactory const& eb_box_factory = warpx.fieldEBFactory(lev);
    amrex::FabArray<amrex::EBCellFlagFab> const& eb_flag = eb_box_factory.getMultiEBCellFlagFab();
    amrex::MultiCutFab const& eb_bnd_cent = eb_box_factory.getBndryCent();
    amrex::MultiCutFab const& eb_bnd_normal = eb_box_factory.getBndryNormal();
    amrex::Array<const amrex::MultiCutFab*,AMREX_SPACEDIM> eb_area_fraction = eb_box_factory.getAreaFrac();

    // get surface integration element
    auto cell_size = WarpX::CellSize(lev);
    amrex::Real const dSx = cell_size[1]*cell_size[2];
    amrex::Real const dSy = cell_size[2]*cell_size[0];
    amrex::Real const dSz = cell_size[0]*cell_size[1];

    // Integral to calculate
    amrex::Gpu::Buffer<amrex::Real> surface_integral({0.0_rt});
    amrex::Real* surface_integral_pointer = surface_integral.data();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // Loop over boxes
    for (amrex::MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box & box = enclosedCells(mfi.nodaltilebox());

        // Skip boxes that are do not intersect with the embedded boundary
        // (i.e. either fully covered or fully regular)
        amrex::FabType fab_type = eb_flag[mfi].getType(box);
        if (fab_type == amrex::FabType::regular) continue;
        if (fab_type == amrex::FabType::covered) continue;

        // Extract data for electric field
        const amrex::Array4<const amrex::Real> & Ex_arr = Ex.array(mfi);
        const amrex::Array4<const amrex::Real> & Ey_arr = Ey.array(mfi);
        const amrex::Array4<const amrex::Real> & Ez_arr = Ez.array(mfi);

        // Extract data for EB
        auto const& eb_flag_arr = eb_flag.array(mfi);
        const amrex::Array4<const amrex::Real> & eb_bnd_normal_arr = eb_bnd_normal.array(mfi);
        const amrex::Array4<const amrex::Real> & eb_bnd_cent_arr = eb_bnd_cent.array(mfi);
        const amrex::Array4<const amrex::Real> & dSx_fraction_arr = eb_area_fraction[0]->array(mfi);
        const amrex::Array4<const amrex::Real> & dSy_fraction_arr = eb_area_fraction[1]->array(mfi);
        const amrex::Array4<const amrex::Real> & dSz_fraction_arr = eb_area_fraction[2]->array(mfi);

        amrex::ParallelFor( box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                // Only cells that are partially covered do contribute to the integral
                if (eb_flag_arr(i,j,k).isRegular() || eb_flag_arr(i,j,k).isCovered()) return;

                // Find nodal point which is outside of the EB
                // (eb_normal points towards the *interior* of the EB)
                int const i_n = (eb_bnd_normal_arr(i,j,k,0) > 0)? i : i+1;
                int const j_n = (eb_bnd_normal_arr(i,j,k,1) > 0)? j : j+1;
                int const k_n = (eb_bnd_normal_arr(i,j,k,2) > 0)? k : k+1;

                // Find cell-centered point which is outside of the EB
                // (eb_normal points towards the *interior* of the EB)
                int i_c = i;
                if ((eb_bnd_normal_arr(i,j,k,0)>0) && (eb_bnd_cent_arr(i,j,k,0)<=0)) i_c -= 1;
                if ((eb_bnd_normal_arr(i,j,k,0)<0) && (eb_bnd_cent_arr(i,j,k,0)>=0)) i_c += 1;
                int j_c = j;
                if ((eb_bnd_normal_arr(i,j,k,1)>0) && (eb_bnd_cent_arr(i,j,k,1)<=0)) j_c -= 1;
                if ((eb_bnd_normal_arr(i,j,k,1)<0) && (eb_bnd_cent_arr(i,j,k,1)>=0)) j_c += 1;
                int k_c = k;
                if ((eb_bnd_normal_arr(i,j,k,2)>0) && (eb_bnd_cent_arr(i,j,k,2)<=0)) k_c -= 1;
                if ((eb_bnd_normal_arr(i,j,k,2)<0) && (eb_bnd_cent_arr(i,j,k,2)>=0)) k_c += 1;

                // Compute contribution to the surface integral $\int dS \cdot E$)
                amrex::Real local_integral_contribution = 0;
                local_integral_contribution += Ex_arr(i_c,j_n,k_n)*dSx*dSx_fraction_arr(i,j,k);
                local_integral_contribution += Ey_arr(i_n,j_c,k_n)*dSy*dSy_fraction_arr(i,j,k);
                local_integral_contribution += Ez_arr(i_n,j_n,k_c)*dSz*dSz_fraction_arr(i,j,k);
                amrex::Gpu::Atomic::AddNoRet( surface_integral_pointer, local_integral_contribution );
        });
    }

    // Reduce across MPI ranks
    surface_integral.copyToHost();
    amrex::Real surface_integral_value = *(surface_integral.hostData());
    amrex::ParallelDescriptor::ReduceRealSum( surface_integral_value );

    // save data
    m_data[0] = PhysConst::ep0 * surface_integral_value;
#endif    
}
// end void ChargeInsideBoundary::ComputeDiags
