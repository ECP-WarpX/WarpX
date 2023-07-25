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

//StationFunctor::StationFunctor(const std::array<const amrex::MultiFab* const, 6> arr_mf_src,
StationFunctor::StationFunctor ( const amrex::MultiFab* const mf_src,
                               amrex::Real zlocation, const int lev,
                               amrex::IntVect crse_ratio, const int ncomp)
    //: ComputeDiagFunctor(ncomp, crse_ratio), m_arr_mf_src(arr_mf_src), m_lev(lev), m_z_location(zlocation)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev), m_z_location(zlocation)
{
}

void
StationFunctor::operator ()(amrex::MultiFab& mf_dst, int dcomp, const int i_buffer) const
{
    // 1. First get slice at given z-location
    auto& warpx = WarpX::GetInstance();
    auto geom = warpx.Geom(m_lev);
    std::unique_ptr< amrex::MultiFab > slice = nullptr;
    const int scomp = 0;
    const int slice_dir = 2;
    bool interpolate = true;
    slice = amrex::get_slice_data(slice_dir, m_z_location, *m_mf_src, geom, scomp,
                                  m_mf_src->nComp(), interpolate);

    // Define MF with dmap of dst mf, with all 6 components from the slice generated from m_mf_src
    const int station_index = static_cast<int> ( ( m_z_location - geom.ProbLo(slice_dir))
                                                 / geom.CellSize(slice_dir) );

    amrex::Box slice_box = m_buffer_box;
    slice_box.setSmall(slice_dir, station_index);
    slice_box.setBig(slice_dir, station_index);

    amrex::BoxArray slice_ba(slice_box);
    slice_ba.maxSize( m_max_box_size);

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
                dst_arr(i,j,k_index,n) = src_arr(i,j,k,n);
            });
    }
    slice = nullptr;
    tmp_slice_ptr = nullptr;
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
}
