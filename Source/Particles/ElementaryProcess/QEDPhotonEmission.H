#ifndef QED_PHOTON_EMISSION_H_
#define QED_PHOTON_EMISSION_H_

#include "WarpXConst.H"
#include "WarpXParticleContainer.H"

#include <BreitWheelerEngineWrapper.H>
#include <QuantumSyncEngineWrapper.H>

struct PhotonEmissionFilterFunc
{
    int m_opt_depth_comp;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const PData& ptd, int i) const noexcept
    {
        const auto opt_depth = ptd.m_runtime_rdata[m_opt_depth_comp][i];
        return (opt_depth < 0.0);
    }
};

struct PhotonEmissionTransformFunc
{
    const QuantumSynchrotronGeneratePhotonAndUpdateMomentum& m_emission_functor;

    bool m_target_has_breit_wheeler = false;
    int m_target_opt_depth_comp = 0;
    BreitWheelerGetOpticalDepth m_target_opt_depth_functor;

    PhotonEmissionTransformFunc(
        const QuantumSynchrotronGeneratePhotonAndUpdateMomentum& emission_functor):
        m_emission_functor{emission_functor}
        {}

    void enable_opt_depth_for_target(int target_opt_depth_comp,
     const BreitWheelerGetOpticalDepth&& opt_depth_functor)
    {
        m_target_has_breit_wheeler = true;
        m_target_opt_depth_comp = target_opt_depth_comp;
        m_target_opt_depth_functor = opt_depth_functor;
    }


    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, SrcData& src, int i_src, int i_dst) const noexcept
    {

        const amrex::ParticleReal me = PhysConst::m_e;
        const amrex::ParticleReal one_over_me = 1./me;

        // Compute electric field amplitude in the particle's frame of
        // reference (particularly important when in boosted frame).
        const amrex::ParticleReal w = src.m_runtime_rdata[PIdx::w][i_src];
        const amrex::ParticleReal ux = src.m_runtime_rdata[PIdx::ux][i_src];
        const amrex::ParticleReal uy = src.m_runtime_rdata[PIdx::uy][i_src];
        const amrex::ParticleReal uz = src.m_runtime_rdata[PIdx::uz][i_src];
        const amrex::ParticleReal ex = src.m_runtime_rdata[PIdx::Ex][i_src];
        const amrex::ParticleReal ey = src.m_runtime_rdata[PIdx::Ey][i_src];
        const amrex::ParticleReal ez = src.m_runtime_rdata[PIdx::Ez][i_src];
        const amrex::ParticleReal bx = src.m_runtime_rdata[PIdx::Bx][i_src];
        const amrex::ParticleReal by = src.m_runtime_rdata[PIdx::By][i_src];
        const amrex::ParticleReal bz = src.m_runtime_rdata[PIdx::Bz][i_src];

        auto px = ux*me;
        auto py = uy*me;
        auto pz = uz*me;

        auto g_w = amrex::ParticleReal{0.0};
        auto g_px = amrex::ParticleReal{0.0};
        auto g_py = amrex::ParticleReal{0.0};
        auto g_pz = amrex::ParticleReal{0.0};

        m_emission_functor.operator()<1>(
            &px, &py, &pz,
            ex, ey, ez,
            bx, by, bz,
            w,
            &g_px, &g_py, &g_pz,
            &g_w);

        src.m_runtime_rdata[PIdx::ux][i_src] = px*one_over_me;
        src.m_runtime_rdata[PIdx::uy][i_src] = py*one_over_me;
        src.m_runtime_rdata[PIdx::uz][i_src] = pz*one_over_me;
        dst.m_runtime_rdata[PIdx::w][i_dst] = g_w;
        dst.m_runtime_rdata[PIdx::ux][i_dst] = g_px*one_over_me;
        dst.m_runtime_rdata[PIdx::uy][i_dst] = g_py*one_over_me;
        dst.m_runtime_rdata[PIdx::uz][i_dst] = g_pz*one_over_me;

        if(m_target_has_breit_wheeler){
            dst.m_runtime_rdata[m_target_opt_depth_comp][i_dst] =
                m_target_opt_depth_functor();
        }


    }
};

#endif //QED_PHOTON_EMISSION_H_
