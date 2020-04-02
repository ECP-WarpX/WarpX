#ifndef QED_SCHWINGER_PROCESS_H_
#define QED_SCHWINGER_PROCESS_H_

#include "Particles/WarpXParticleContainer.H"

struct SchwingerFilterFunc
{
    const amrex::Real m_y_size;
    int m_threshold_poisson_gaussian;
    const amrex::Real m_dVdt;

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

        return getSchwingerProductionNumber( m_dVdt,
                    arrEx(i,j,k),arrEy(i,j,k),arrEz(i,j,k),
                    arrBx(i,j,k),arrBy(i,j,k),arrBz(i,j,k),
                    m_threshold_poisson_gaussian);
    }
};

struct SchwingerTransformFunc
{

    const amrex::Real m_y_size;
    const int m_weight_index;

    template <typename DstData1, typename DstData2>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData1& dst1, DstData2& dst2, int i_dst1, int i_dst2,
                               int N, amrex::Real total_weight) const noexcept
    {
        for (int n = 0; n < N; ++n){
#if (AMREX_SPACEDIM == 3)
            dst1.m_rdata[m_weight_index][i_dst1+n] = total_weight/N;
            dst2.m_rdata[m_weight_index][i_dst2+n] = total_weight/N;
#elif (AMREX_SPACEDIM == 2)
            dst1.m_rdata[m_weight_index][i_dst1+n] = total_weight/N/m_y_size;
            dst2.m_rdata[m_weight_index][i_dst2+n] = total_weight/N/m_y_size;
#endif
            }
    }
};

#endif
