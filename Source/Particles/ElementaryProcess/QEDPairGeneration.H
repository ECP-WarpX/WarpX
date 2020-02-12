#ifndef QED_PAIR_GENERATION_H_
#define QED_PAIR_GENERATION_H_

#include "WarpXConst.H"
#include "WarpXParticleContainer.H"

#include <BreitWheelerEngineWrapper.H>
#include <QuantumSyncEngineWrapper.H>


struct PairGenerationFilterFunc
{
    int m_opt_depth_runtime_comp;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const PData& ptd, int i) const noexcept
    {
        const auto opt_depth =
            ptd.m_runtime_rdata[m_opt_depth_runtime_comp][i];
        return (opt_depth < 0.0);
    }
};

struct PairGenerationTransformFunc
{

    const BreitWheelerGeneratePairs& m_generate_functor;

    bool m_target1_has_quantum_sync = false;
    int m_target1_opt_depth_runtime_comp = 0;
    bool m_target2_has_quantum_sync = false;
    int m_target2_opt_depth_runtime_comp = 0;
    QuantumSynchrotronGetOpticalDepth m_target_opt_depth_functor;

    PairGenerationTransformFunc(
        const BreitWheelerGeneratePairs& generate_functor):
        m_generate_functor{generate_functor}
        {}

    void enable_opt_depth_for_targets(const bool enable1, const bool enable2,
     const int target1_opt_depth_runtime_comp, const int target2_opt_depth_runtime_comp,
     const QuantumSynchrotronGetOpticalDepth&& opt_depth_functor)
    {
        if(enable1){
            m_target1_has_quantum_sync = true;
            m_target1_opt_depth_runtime_comp = target1_opt_depth_runtime_comp;
        }

        if(enable2){
            m_target2_has_quantum_sync = true;
            m_target2_opt_depth_runtime_comp = target2_opt_depth_runtime_comp;
        }

        if(enable1 || enable2){
            m_target_opt_depth_functor = opt_depth_functor;
        }

    }

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst1, DstData& dst2, SrcData& src,
        int i_src, int i_dst1, int i_dst2) const noexcept
    {

        const amrex::ParticleReal me = PhysConst::m_e;
        const amrex::ParticleReal one_over_me = 1./me;

        // Compute electric field amplitude in the particle's frame of
        // reference (particularly important when in boosted frame).
        const amrex::ParticleReal w = src.m_rdata[PIdx::w][i_src];
        const amrex::ParticleReal ux = src.m_rdata[PIdx::ux][i_src];
        const amrex::ParticleReal uy = src.m_rdata[PIdx::uy][i_src];
        const amrex::ParticleReal uz = src.m_rdata[PIdx::uz][i_src];
        const amrex::ParticleReal ex = src.m_rdata[PIdx::Ex][i_src];
        const amrex::ParticleReal ey = src.m_rdata[PIdx::Ey][i_src];
        const amrex::ParticleReal ez = src.m_rdata[PIdx::Ez][i_src];
        const amrex::ParticleReal bx = src.m_rdata[PIdx::Bx][i_src];
        const amrex::ParticleReal by = src.m_rdata[PIdx::By][i_src];
        const amrex::ParticleReal bz = src.m_rdata[PIdx::Bz][i_src];

        const auto px = ux*me;
        const auto py = uy*me;
        const auto pz = uz*me;

        auto e_w = amrex::ParticleReal{0.0};
        auto p_w = amrex::ParticleReal{0.0};
        auto e_px = amrex::ParticleReal{0.0};
        auto e_py = amrex::ParticleReal{0.0};
        auto e_pz = amrex::ParticleReal{0.0};
        auto p_px = amrex::ParticleReal{0.0};
        auto p_py = amrex::ParticleReal{0.0};
        auto p_pz = amrex::ParticleReal{0.0};

        m_generate_functor.operator()<1>(
            px, py, pz,
            ex, ey, ez,
            bx, by, bz,
            w,
            &e_px, &e_py, &e_pz,
            &p_px, &p_py, &p_pz,
            &e_w, &p_w);

        src.m_rdata[PIdx::ux][i_src] = px*one_over_me;
        src.m_rdata[PIdx::uy][i_src] = py*one_over_me;
        src.m_rdata[PIdx::uz][i_src] = pz*one_over_me;

        dst1.m_rdata[PIdx::w][i_dst1] = e_w;
        dst1.m_rdata[PIdx::ux][i_dst1] = e_px*one_over_me;
        dst1.m_rdata[PIdx::uy][i_dst1] = e_py*one_over_me;
        dst1.m_rdata[PIdx::uz][i_dst1] = e_pz*one_over_me;

        dst2.m_rdata[PIdx::w][i_dst2] = p_w;
        dst2.m_rdata[PIdx::ux][i_dst2] = p_px*one_over_me;
        dst2.m_rdata[PIdx::uy][i_dst2] = p_py*one_over_me;
        dst2.m_rdata[PIdx::uz][i_dst2] = p_pz*one_over_me;

        if(m_target1_has_quantum_sync){
            dst1.m_runtime_rdata[m_target1_opt_depth_runtime_comp][i_dst1] =
                m_target_opt_depth_functor();
        }

        if(m_target2_has_quantum_sync){
            dst2.m_runtime_rdata[m_target2_opt_depth_runtime_comp][i_dst2] =
                m_target_opt_depth_functor();
        }

    }

};

#endif //QED_PAIR_GENERATION_H_
