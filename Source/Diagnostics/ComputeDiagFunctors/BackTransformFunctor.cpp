/* Copyright 2021 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackTransformFunctor.H"

#include "Diagnostics/ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/utils/Communication.H>

#include <AMReX_Array4.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include <cmath>
#include <map>
#include <memory>

using namespace amrex;

BackTransformFunctor::BackTransformFunctor (amrex::MultiFab const * mf_src, int lev,
                                            const int ncomp, const int num_buffers,
                                            amrex::Vector< std::string > varnames,
                                            const amrex::IntVect crse_ratio)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev), m_num_buffers(num_buffers), m_varnames(varnames)
{
    InitData();
}

void
BackTransformFunctor::operator ()(amrex::MultiFab& mf_dst, int /*dcomp*/, const int i_buffer) const
{
    // Perform back-transformation only if z slice is within the domain stored as 0/1
    // in m_perform_backtransform[i_buffer]
    if ( m_perform_backtransform[i_buffer] == 1) {
        auto& warpx = WarpX::GetInstance();
        auto geom = warpx.Geom(m_lev);
        amrex::Real gamma_boost = warpx.gamma_boost;
        int moving_window_dir = warpx.moving_window_dir;
        amrex::Real beta_boost = std::sqrt( 1._rt - 1._rt/( gamma_boost * gamma_boost) );
        const bool interpolate = true;
        std::unique_ptr< amrex::MultiFab > slice = nullptr;
        int scomp = 0;
        // Generate slice of the cell-centered multifab containing boosted-frame field-data
        // at current z-boost location for the ith buffer
        slice = amrex::get_slice_data(moving_window_dir,
                                     m_current_z_boost[i_buffer],
                                     *m_mf_src,
                                     geom,
                                     scomp,
                                     m_mf_src->nComp(),
                                     interpolate);
        // Perform in-place Lorentz-transform of all the fields stored in the slice.
        LorentzTransformZ( *slice, gamma_boost, beta_boost);

        // Create a 2D box for the slice in the boosted frame
        amrex::Real dx = geom.CellSize(moving_window_dir);
        // index corresponding to z_boost location in the boost-frame
        int i_boost = static_cast<int> ( ( m_current_z_boost[i_buffer]
                                            - geom.ProbLo(moving_window_dir) ) / dx );
        // z-Slice at i_boost with x,y indices same as buffer_box
        amrex::Box slice_box = m_buffer_box[i_buffer];
        slice_box.setSmall(moving_window_dir, i_boost);
        slice_box.setBig(moving_window_dir, i_boost);

        // Make it a BoxArray
        amrex::BoxArray slice_ba(slice_box);
        slice_ba.maxSize( m_max_box_size );
        // Define MultiFab with the distribution map of the destination multifab and
        // containing all ten components that were in the slice generated from m_mf_src.
        std::unique_ptr< amrex::MultiFab > tmp_slice_ptr = nullptr;
        tmp_slice_ptr = std::make_unique<MultiFab> ( slice_ba, mf_dst.DistributionMap(),
                                                     slice->nComp(), 0 );
        tmp_slice_ptr->setVal(0.0);
        // Parallel copy the lab-frame data from "slice" MultiFab with
        // ncomp=10 and boosted-frame dmap to "tmp_slice_ptr" MultiFab with
        // ncomp=10 and dmap of the destination Multifab, which will store the final data
        ablastr::utils::communication::ParallelCopy(*tmp_slice_ptr, *slice, 0, 0, slice->nComp(),
                                                    IntVect(AMREX_D_DECL(0, 0, 0)),
                                                    IntVect(AMREX_D_DECL(0, 0, 0)),
                                                    WarpX::do_single_precision_comms);
        // Now we will cherry pick only the user-defined fields from
        // tmp_slice_ptr to dst_mf
        const int k_lab = m_k_index_zlab[i_buffer];
        const int ncomp_dst = mf_dst.nComp();
        amrex::MultiFab& tmp = *tmp_slice_ptr;
#ifdef AMREX_USE_GPU
        Gpu::DeviceVector<int> d_map_varnames(m_map_varnames.size());
        Gpu::copyAsync(Gpu::hostToDevice,
                       m_map_varnames.begin(), m_map_varnames.end(),
                       d_map_varnames.begin());
        Gpu::synchronize();
        int const* field_map_ptr = d_map_varnames.dataPtr();
#else
        int const* field_map_ptr = m_map_varnames.dataPtr();
#endif
        for (amrex::MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            amrex::Array4<amrex::Real> src_arr = tmp[mfi].array();
            amrex::Array4<amrex::Real> dst_arr = mf_dst[mfi].array();
            amrex::ParallelFor( tbx, ncomp_dst,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
                {
                    const int icomp = field_map_ptr[n];
#if defined(WARPX_DIM_3D)
                    dst_arr(i, j, k_lab, n) = src_arr(i, j, k, icomp);
#else
                    dst_arr(i, k_lab, k, n) = src_arr(i, j, k, icomp);
#endif
                } );
        }

        // Reset the temporary MultiFabs generated
        slice = nullptr;
        tmp_slice_ptr = nullptr;
    }

}

void
BackTransformFunctor::PrepareFunctorData (int i_buffer,
                          bool z_slice_in_domain, amrex::Real current_z_boost,
                          amrex::Box buffer_box, const int k_index_zlab,
                          const int max_box_size, const int snapshot_full)
{
    m_buffer_box[i_buffer] = buffer_box;
    m_current_z_boost[i_buffer] = current_z_boost;
    m_k_index_zlab[i_buffer] = k_index_zlab;
    m_perform_backtransform[i_buffer] = 0;
    if (z_slice_in_domain == true and snapshot_full == 0) m_perform_backtransform[i_buffer] = 1;
    m_max_box_size = max_box_size;
}

void
BackTransformFunctor::InitData ()
{

    m_buffer_box.resize( m_num_buffers );
    m_current_z_boost.resize( m_num_buffers );
    m_perform_backtransform.resize( m_num_buffers );
    m_k_index_zlab.resize( m_num_buffers );
    m_map_varnames.resize( m_varnames.size() );

    std::map<std::string, int> m_possible_fields_to_dump = {
        {"Ex", 0},
        {"Ey", 1},
        {"Ez", 2},
        {"Bx", 3},
        {"By", 4},
        {"Bz", 5},
        {"jx", 6},
        {"jy", 7},
        {"jz", 8},
        {"rho", 9}
    };

    for (int i = 0; i < m_varnames.size(); ++i)
    {
        m_map_varnames[i] = m_possible_fields_to_dump[ m_varnames[i] ] ;
    }

}

void
BackTransformFunctor::LorentzTransformZ (amrex::MultiFab& data, amrex::Real gamma_boost,
                                         amrex::Real beta_boost) const
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(data, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& tbx = mfi.tilebox();
        amrex::Array4< amrex::Real > arr = data[mfi].array();
        amrex::Real clight = PhysConst::c;
        amrex::Real inv_clight = 1.0_rt/clight;
        // arr(x,y,z,comp) has ten-components namely,
        // Ex Ey Ez Bx By Bz jx jy jz rho in that order.
        amrex::ParallelFor( tbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Back-transform the transverse electric and magnetic fields.
                // Note that the z-components, Ez, Bz, are not changed by the transform.
                amrex::Real e_lab, b_lab, j_lab, rho_lab;
                // Transform Ex_boost (ncomp=0) & By_boost (ncomp=4) to lab-frame
                e_lab = gamma_boost * ( arr(i, j, k, 0)
                                        + beta_boost * clight * arr(i, j, k, 4) );
                b_lab = gamma_boost * ( arr(i, j, k, 4)
                                        + beta_boost * inv_clight * arr(i, j, k, 0) );
                // Store lab-frame data in-place
                arr(i, j, k, 0) = e_lab;
                arr(i, j, k, 4) = b_lab;

                // Transform Ey_boost (ncomp=1) & Bx_boost (ncomp=3) to lab-frame
                e_lab = gamma_boost * ( arr(i, j, k, 1)
                                        - beta_boost * clight * arr(i, j, k, 3) );
                b_lab = gamma_boost * ( arr(i, j, k, 3)
                                        - beta_boost * inv_clight * arr(i, j, k, 1) );
                // Store lab-frame data in-place
                arr(i, j, k, 1) = e_lab;
                arr(i, j, k, 3) = b_lab;

                // Transform charge density (ncomp=9)
                // and z-component of current density (ncomp=8)
                j_lab = gamma_boost * ( arr(i, j, k, 8)
                                        + beta_boost * clight * arr(i, j, k, 9) );
                rho_lab = gamma_boost * ( arr(i, j, k, 9)
                                          + beta_boost * inv_clight * arr(i, j, k, 8) );
                // Store lab-frame jz and rho in-place
                arr(i, j, k, 8) = j_lab;
                arr(i, j, k, 9) = rho_lab;
            }
        );
    }

}
