#ifndef QED_PAIR_GENERATION_H_
#define QED_PAIR_GENERATION_H_

#include "WarpXConst.H"
#include "WarpXParticleContainer.H"

struct PhotonEmissionFilterFunc
{
    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const PData& ptd, int i) const noexcept
    {

        return false;
    }
};

struct PhotonEmissionTransformFunc
{
    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, SrcData& src, int i_src, int i_dst) const noexcept
    {

    }
};

#endif //QED_PAIR_GENERATION_H_
