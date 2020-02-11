#ifndef QED_PHOTON_EMISSION_H_
#define QED_PHOTON_EMISSION_H_

#include "WarpXConst.H"
#include "WarpXParticleContainer.H"

#include <QuantumSyncEngineWrapper.H>

struct PhotonEmissionFilterFunc
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

struct PhotonEmissionTransformFunc
{
    const QuantumSynchrotronGeneratePhotonAndUpdateMomentum m_emission_functor;

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

        if(sqrt(g_px*g_px + g_py*g_py + g_pz*g_pz) < 2.0 * me*PhysConst::c){
            //p_product_phot.id() = -1;
        }

        /*
                    r_source_attribs[PIdx::ux][is] = px*one_over_me;
                    r_source_attribs[PIdx::uy][is] = py*one_over_me;
                    r_source_attribs[PIdx::uz][is] = pz*one_over_me;
                    r_phot_attribs[PIdx::w][ip] = g_w;
                    r_phot_attribs[PIdx::ux][ip]= g_px*one_over_me;
                    r_phot_attribs[PIdx::uy][ip] = g_py*one_over_me;
                    r_phot_attribs[PIdx::uz][ip] = g_pz*one_over_me;
        */

    }
};

#endif //QED_PHOTON_EMISSION_H_
