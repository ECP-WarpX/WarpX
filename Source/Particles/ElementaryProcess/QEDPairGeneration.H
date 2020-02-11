#ifndef QED_PAIR_GENERATION_H_
#define QED_PAIR_GENERATION_H_

#include "WarpXConst.H"
#include "WarpXParticleContainer.H"

#include <BreitWheelerEngineWrapper.H>


struct PairGenerationFilterFunc
{
    int m_opt_depth_comp;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const PData& ptd, int i) const noexcept
    {
        const auto opt_depth = ptd.m_runtime_data[m_opt_depth_comp][i];
        return (opt_depth < 0.0);
    }
};

struct PairGenerationTransformFunc
{

    const BreitWheelerGeneratePairs m_generate_functor;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, SrcData& src, int i_src, int i_dst) const noexcept
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
            /*
                                r_ele_attribs[PIdx::w][ip] = e_w;
                    r_ele_attribs[PIdx::ux][ip]= e_px*one_over_me;
                    r_ele_attribs[PIdx::uy][ip] = e_py*one_over_me;
                    r_ele_attribs[PIdx::uz][ip] = e_pz*one_over_me;
                    r_pos_attribs[PIdx::w][ip] = p_w;
                    r_pos_attribs[PIdx::ux][ip]= p_px*one_over_me;
                    r_pos_attribs[PIdx::uy][ip] = p_py*one_over_me;
                    r_pos_attribs[PIdx::uz][ip] = p_pz*one_over_me;

                    p_source.id() = -1; //Sets ID to zero */
    }

};

#endif //QED_PAIR_GENERATION_H_
