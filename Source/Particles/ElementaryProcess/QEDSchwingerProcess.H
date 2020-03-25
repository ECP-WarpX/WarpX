#ifndef QED_SCHWINGER_PROCESS_H_
#define QED_SCHWINGER_PROCESS_H_

struct SchwingerFilterFunc
{
    const amrex::Real m_y_size;
    int m_threshold_poisson_gaussian;

    template <typename FABs>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (const FABs *src_FABs, int i, int j, int k) const noexcept
    {
        const auto& arrEx = src_FABs[0]->array();
        const auto& arrEy = src_FABs[1]->array();
        const auto& arrEz = src_FABs[2]->array();
        const auto& arrBx = src_FABs[3]->array();
        const auto& arrBy = src_FABs[4]->array();
        const auto& arrBz = src_FABs[5]->array();

    // So far Schwinger process does not work with mesh refinement
    const int level_0 = 0;

    // get cell size multiplied by time step
    auto & warpx = WarpX::GetInstance();
    amrex::Geometry const & geom = warpx.Geom(level_0);
    auto domain_box = geom.Domain();
#if (AMREX_SPACEDIM == 2)
    auto dVdt = geom.CellSize(0) * geom.CellSize(1) * m_y_size
               * warpx.getdt(level_0);
#elif (AMREX_SPACEDIM == 3)
    auto dVdt = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2)
               * warpx.getdt(level_0);
#endif

        return getSchwingerProductionNumber( dVdt,
                    arrEx(i,j,k),arrEy(i,j,k),arrEz(i,j,k),
                    arrBx(i,j,k),arrBy(i,j,k),arrBz(i,j,k),
                    m_threshold_poisson_gaussian);
    }
};

struct SchwingerTransformFunc
{
    template <typename DstData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, int i_dst, int N, amrex::Real total_weight) const noexcept
    {
    }
};

#endif
