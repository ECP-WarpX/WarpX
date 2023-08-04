#include "StationFunctor.H"

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

//StationFunctor::StationFunctor ( const amrex::MultiFab* const mf_src,
StationFunctor::StationFunctor(const std::array<const amrex::MultiFab* const, 6> arr_mf_src,
                               amrex::Real zlocation, const int lev,
                               amrex::IntVect crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_arr_mf_src(arr_mf_src), m_lev(lev), m_z_location(zlocation)
//    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev), m_z_location(zlocation)
{
}

void
StationFunctor::operator ()(amrex::MultiFab& mf_dst, const int dcomp, const int i_buffer) const
{
    // 1. First get slice at given z-location
    if (m_slice_in_domain == 1) {
        auto& warpx = WarpX::GetInstance();
        auto geom = warpx.Geom(m_lev);
//        std::unique_ptr< amrex::MultiFab > slice = nullptr;
        const int scomp = 0;
        const int slice_dir = WARPX_ZINDEX;
        bool interpolate = true;

        amrex::Vector<int> slice_to_full_ba_map;
        std::unique_ptr<amrex::MultiFab> slice = AllocateSlice(slice_to_full_ba_map);
        slice->setVal(0.);

        for (amrex::MFIter mfi(*slice, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            int slice_gid = mfi.index();
            int full_gid = slice_to_full_ba_map[slice_gid];
            const amrex::Array4<amrex::Real> slice_arr = (*slice)[mfi].array();
            const amrex::Array4<amrex::Real const> src_Ex_arr = (*m_arr_mf_src[0])[full_gid].array();
            const amrex::Array4<amrex::Real const> src_Ey_arr = (*m_arr_mf_src[1])[full_gid].array();
            const amrex::Array4<amrex::Real const> src_Ez_arr = (*m_arr_mf_src[2])[full_gid].array();
            const amrex::Array4<amrex::Real const> src_Bx_arr = (*m_arr_mf_src[3])[full_gid].array();
            const amrex::Array4<amrex::Real const> src_By_arr = (*m_arr_mf_src[4])[full_gid].array();
            const amrex::Array4<amrex::Real const> src_Bz_arr = (*m_arr_mf_src[5])[full_gid].array();

            const amrex::Box& tbx = mfi.tilebox();
            amrex::ParallelFor( tbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    slice_arr(i,j,k,0) = src_Ex_arr(i,j,k);
                    slice_arr(i,j,k,1) = src_Ey_arr(i,j,k);
                    slice_arr(i,j,k,2) = src_Ez_arr(i,j,k);
                    slice_arr(i,j,k,3) = src_Bx_arr(i,j,k);
                    slice_arr(i,j,k,4) = src_By_arr(i,j,k);
                    slice_arr(i,j,k,5) = src_Bz_arr(i,j,k);
                }
            );
        }

        // Define MF with dmap of dst mf, with all 6 components from the slice generated from m_mf_src
        const int station_index = static_cast<int> ( ( m_z_location - geom.ProbLo(slice_dir))
                                                     / geom.CellSize(slice_dir) );
        amrex::Box slice_box = amrex::surroundingNodes(m_buffer_box);
        slice_box.setSmall(slice_dir, station_index);
        slice_box.setBig(slice_dir, station_index);

        amrex::BoxArray slice_ba(slice_box);
        slice_ba.maxSize( 256);
        std::unique_ptr< amrex::MultiFab > tmp_slice_ptr = nullptr;
        const int nghost = 1;
        tmp_slice_ptr = std::make_unique< amrex::MultiFab > (slice_ba, mf_dst.DistributionMap(),
                                                             slice->nComp(), nghost);
        tmp_slice_ptr->setVal(0.);

        // Parallel copy slice to tmp_slice MF
        const int dcomp_tmp = 0;
        amrex::IntVect src_ngrow = amrex::IntVect( AMREX_D_DECL(1,1,1) );
        amrex::IntVect dst_ngrow = amrex::IntVect( AMREX_D_DECL(1,1,1) );
        ablastr::utils::communication::ParallelCopy(*tmp_slice_ptr, *slice, scomp, dcomp_tmp, slice->nComp(),
                                                    slice->nGrowVect(), tmp_slice_ptr->nGrowVect(),
                                                    WarpX::do_single_precision_comms);

        // MFIter to copy tmp_slice at k-index location in time-based mf
        const int k_index = m_k_index;
        const int ncomp_dst = mf_dst.nComp();
        amrex::MultiFab& tmp = *tmp_slice_ptr;
        for (amrex::MFIter mfi(tmp, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box & tbx = mfi.tilebox();
            const amrex::Array4<amrex::Real> src_arr = tmp[mfi].array();
            const amrex::Array4<amrex::Real> dst_arr = mf_dst[mfi].array();
            amrex::ParallelFor( tbx, ncomp_dst,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
#if defined(WARPX_DIM_3D)
                    dst_arr(i,j,k_index,n) = src_arr(i,j,k,n);
#elif defined(WARPX_DIM_XZ)
                    dst_arr(i,k_index,k,n) = src_arr(i,j,k,n);
#endif
                });
        }
        slice = nullptr;
        tmp_slice_ptr = nullptr;
    }
}

void
StationFunctor::PrepareFunctorData (int i_station, bool slice_in_domain, amrex::Real zlocation,
                                    amrex::Box buffer_box, int k_index, const int max_box_size,
                                    const int buffer_full)
{
    amrex::ignore_unused(i_station, slice_in_domain, zlocation, buffer_full);
    m_buffer_box = buffer_box;
    m_k_index = k_index;
    m_max_box_size = max_box_size;
    if (slice_in_domain == true) {
        m_slice_in_domain = 1;
    } else {
        m_slice_in_domain = 0;
    }
}

std::unique_ptr<amrex::MultiFab>
StationFunctor::AllocateSlice (amrex::Vector<int>& slice_to_full_ba_map) const
{
    auto& warpx = WarpX::GetInstance();
    amrex::Geometry geom = warpx.Geom(m_lev);
    // get slice and convert to index space
    amrex::RealBox real_slice = geom.ProbDomain();
    real_slice.setLo(WARPX_ZINDEX, m_z_location);
    real_slice.setHi(WARPX_ZINDEX, m_z_location);
    amrex::IntVect slice_lo = amrex::IntVect( AMREX_D_DECL(
                              static_cast<int>(std::floor((real_slice.lo(0) - geom.ProbLo(0))/geom.CellSize(0))),
                              static_cast<int>(std::floor((real_slice.lo(1) - geom.ProbLo(1))/geom.CellSize(1))),
                              static_cast<int>(std::floor((real_slice.lo(2) - geom.ProbLo(2))/geom.CellSize(2)))
                              ));
    amrex::IntVect slice_hi = amrex::IntVect( AMREX_D_DECL(
                              static_cast<int>(std::floor((real_slice.hi(0) - geom.ProbLo(0))/geom.CellSize(0))),
                              static_cast<int>(std::floor((real_slice.hi(1) - geom.ProbLo(1))/geom.CellSize(1))),
                              static_cast<int>(std::floor((real_slice.hi(2) - geom.ProbLo(2))/geom.CellSize(2)))
                             ));
    amrex::Box slice_box = amrex::surroundingNodes(amrex::Box(slice_lo, slice_hi));

    // Define nodal multifab that stores the slice
    amrex::BoxArray const& ba = amrex::convert( m_arr_mf_src[0]->boxArray(), amrex::IntVect::TheNodeVector());
    const amrex::DistributionMapping& dm = m_arr_mf_src[0]->DistributionMap();
    std::vector< std::pair<int, amrex::Box> > isects;
    ba.intersections(slice_box, isects, false, 0);
    amrex::Vector<amrex::Box> boxes;
    amrex::Vector<int> procs;
    for (auto const& is : isects) {
        procs.push_back(dm[is.first]);
        boxes.push_back(is.second);
        slice_to_full_ba_map.push_back(is.first);
    }
    amrex::BoxArray slice_ba(&boxes[0], static_cast<int>(boxes.size()));
    amrex::DistributionMapping slice_dmap(std::move(procs));
    const int nghost = 1;
    std::unique_ptr<amrex::MultiFab> slice = std::make_unique<amrex::MultiFab>(slice_ba, slice_dmap, nComp(), nghost,
                                                                     amrex::MFInfo(), amrex::FArrayBoxFactory());
    return slice;
}
