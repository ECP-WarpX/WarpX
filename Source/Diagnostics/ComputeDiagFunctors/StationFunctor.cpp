#include "StationFunctor.H"

#include "WarpX.H"

#include <ablastr/utils/Communication.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <algorithm>

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real interp_to_slice (int i, int j,
                             amrex::Array4<amrex::Real const> const& src,
                             amrex::IndexType type, amrex::Real z)
{
    using amrex::Real;
    if (type.cellCentered()) { z -= Real(0.5); }
    int const iz = static_cast<int>(std::floor(z));
    Real const w = z - Real(iz);
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(i,j);
    amrex::Abort("xxxxx todo");
    return 0;
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(j);
    if (type.nodeCentered(0)) {
        return src(i,iz,0)*(Real(1.0)-w) + src(i,iz+1,0)*w;
    } else {
        return Real(0.5)*(src(i-1,iz,0)*(Real(1.0)-w) + src(i-1,iz+1,0)*w +
                          src(i  ,iz,0)*(Real(1.0)-w) + src(i  ,iz+1,0)*w);
    }
#else
    amrex::Abort("xxxxx todo");
    return 0;
#endif
}

}

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
    if (! m_slice_in_domain) { return; }

    // 1. First get slice at given z-location
    auto const& warpx = WarpX::GetInstance();
    auto const& geom = warpx.Geom(m_lev);
    amrex::Box const& domain = geom.Domain();
    const int slice_dir = WARPX_ZINDEX;
    // cell index for the slice
    int slice_k = static_cast<int>(std::floor((m_z_location - geom.ProbLo(slice_dir))
                                              * geom.InvCellSize(slice_dir)))
        + domain.smallEnd(slice_dir);
    slice_k = std::clamp(slice_k, domain.smallEnd(WARPX_ZINDEX), domain.bigEnd(WARPX_ZINDEX));
    const amrex::Real slice_z = (m_z_location - geom.ProbLo(slice_dir)) *
        geom.InvCellSize(slice_dir) + amrex::Real(domain.smallEnd(slice_dir));

    amrex::Vector<int> slice_to_full_ba_map;
    std::unique_ptr<amrex::MultiFab> slice = AllocateSlice(slice_to_full_ba_map, slice_k);

    for (auto& smf : m_arr_mf_src) {
        AMREX_ALWAYS_ASSERT(smf->nGrowVect().allGE(amrex::IntVect(1)));
        const_cast<amrex::MultiFab*>(smf)->FillBoundary
            (0,1,amrex::IntVect(1),geom.periodicity());
    }

    amrex::GpuArray<amrex::IndexType,6> src_itype
        {m_arr_mf_src[0]->ixType(),
         m_arr_mf_src[1]->ixType(),
         m_arr_mf_src[2]->ixType(),
         m_arr_mf_src[3]->ixType(),
         m_arr_mf_src[4]->ixType(),
         m_arr_mf_src[5]->ixType()};

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*slice, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        int full_gid = slice_to_full_ba_map[mfi.LocalIndex()];
        const amrex::Array4<amrex::Real> slice_arr = (*slice)[mfi].array();
        amrex::GpuArray<amrex::Array4<amrex::Real const>,6> src_arr
            {(*m_arr_mf_src[0])[full_gid].const_array(),
             (*m_arr_mf_src[1])[full_gid].const_array(),
             (*m_arr_mf_src[2])[full_gid].const_array(),
             (*m_arr_mf_src[3])[full_gid].const_array(),
             (*m_arr_mf_src[4])[full_gid].const_array(),
             (*m_arr_mf_src[5])[full_gid].const_array()};

        const amrex::Box& tbx = mfi.tilebox();
        amrex::ParallelFor( tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            for (int n = 0; n < 6; ++n) {
                slice_arr(i,j,k,n) = interp_to_slice(i,j,src_arr[n],src_itype[n],slice_z);
            }
        });
    }

    mf_dst.ParallelCopy(*slice, 0, dcomp, nComp());
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
StationFunctor::AllocateSlice (amrex::Vector<int>& slice_to_full_ba_map, int slice_k) const
{
    auto const& warpx = WarpX::GetInstance();
    amrex::Geometry const& geom = warpx.Geom(m_lev);
    amrex::Box slice_box = geom.Domain();
    slice_box.setSmall(WARPX_ZINDEX, slice_k).setBig(WARPX_ZINDEX, slice_k);

    amrex::BoxArray const& ba = amrex::convert(m_arr_mf_src[0]->boxArray(),
                                               amrex::IntVect::TheCellVector());
    const amrex::DistributionMapping& dm = m_arr_mf_src[0]->DistributionMap();
    std::vector< std::pair<int, amrex::Box> > isects;
    ba.intersections(slice_box, isects);
    amrex::BoxList boxes;
    amrex::Vector<int> procs;
    for (auto const& is : isects) {
        procs.push_back(dm[is.first]);
        auto b = is.second;
        b.setSmall(WARPX_ZINDEX, m_k_index).setBig(WARPX_ZINDEX, m_k_index);
        boxes.push_back(b);
        if (dm[is.first] == amrex::ParallelDescriptor::MyProc()) {
            slice_to_full_ba_map.push_back(is.first);
        }
    }
    amrex::BoxArray slice_ba(std::move(boxes));
    amrex::IntVect typ(1);
    typ[WARPX_ZINDEX] = 0; // nodal except in z-direction
    slice_ba.convert(typ);
    amrex::DistributionMapping slice_dmap(std::move(procs));
    const int nghost = 0;
    return std::make_unique<amrex::MultiFab>(slice_ba, slice_dmap, nComp(), nghost);
}
