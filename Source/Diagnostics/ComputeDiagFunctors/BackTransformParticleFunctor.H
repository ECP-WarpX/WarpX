#ifndef WARPX_BACKTRANSFORMPARTICLEFUNCTOR_H_
#define WARPX_BACKTRANSFORMPARTICLEFUNCTOR_H_

#include "ComputeParticleDiagFunctor.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <AMReX.H>
#include <AMReX_AmrParticles.H>





struct SelectParticles
{
    using TmpParticles = WarpXParticleContainer::TmpParticles;

    SelectParticles( const WarpXParIter& a_pti, TmpParticles& tmp_particle_data,
                     amrex::Real current_z_boost, amrex::Real old_z_boost,
                     int a_offset = 0);

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
        amrex::ignore_unused(src);
        using namespace amrex;
        amrex::ParticleReal xp, yp, zp;
        m_get_position(i, xp, yp, zp);
        int Flag = 0;
        if ( ( (zp >= m_current_z_boost) && (zpold[i] <= m_old_z_boost) ) ||
             ( (zp <= m_current_z_boost) && (zpold[i] >= m_old_z_boost) ))
        {    Flag = 1;
        }
        return Flag;
    }

    GetParticlePosition m_get_position;
    amrex::Real m_current_z_boost;
    amrex::Real m_old_z_boost;
    amrex::ParticleReal* AMREX_RESTRICT zpold = nullptr;
};

struct LorentzTransformParticles
{

    using TmpParticles = WarpXParticleContainer::TmpParticles;

    LorentzTransformParticles ( const WarpXParIter& a_pti, TmpParticles& tmp_particle_data,
                                amrex::Real t_boost, amrex::Real dt,
                                amrex::Real t_lab, int a_offset = 0);

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator () (const DstData& dst, const SrcData& src, int i_src, int i_dst) const noexcept
    {
        amrex::ignore_unused(src);
        using namespace amrex::literals;
        // get current src position
        amrex::ParticleReal xpnew, ypnew, zpnew;
        m_get_position(i_src, xpnew, ypnew, zpnew);
        const amrex::Real gamma_new_p = std::sqrt(1.0_rt + m_inv_c2*
                                        ( m_uxpnew[i_src] * m_uxpnew[i_src]
                                        + m_uypnew[i_src] * m_uypnew[i_src]
                                        + m_uzpnew[i_src] * m_uzpnew[i_src]));
        const amrex::Real gamma_old_p = std::sqrt(1.0_rt + m_inv_c2*
                                        ( m_uxpold[i_src] * m_uxpold[i_src]
                                        + m_uypold[i_src] * m_uypold[i_src]
                                        + m_uzpold[i_src] * m_uzpold[i_src]));
        const amrex::Real t_new_p = m_gammaboost * m_t_boost - m_uzfrm * zpnew * m_inv_c2;
        const amrex::Real z_new_p = m_gammaboost* ( zpnew + m_betaboost * m_Phys_c * m_t_boost);
        const amrex::Real uz_new_p = m_gammaboost * m_uzpnew[i_src] - gamma_new_p * m_uzfrm;
        const amrex::Real t_old_p = m_gammaboost * (m_t_boost - m_dt)
                                    - m_uzfrm * m_zpold[i_src] * m_inv_c2;
        const amrex::Real z_old_p = m_gammaboost * ( m_zpold[i_src] + m_betaboost
                                                     * m_Phys_c * (m_t_boost - m_dt ) );
        const amrex::Real uz_old_p = m_gammaboost * m_uzpold[i_src] - gamma_old_p * m_uzfrm;
        // interpolate in time to t_lab
        const amrex::Real weight_old = (t_new_p - m_t_lab)
                                     / (t_new_p - t_old_p);
        const amrex::Real weight_new = (m_t_lab - t_old_p)
                                     / (t_new_p - t_old_p);
        // weighted sum of old and new values
        const amrex::ParticleReal xp = m_xpold[i_src] * weight_old + xpnew * weight_new;
        const amrex::ParticleReal yp = m_ypold[i_src] * weight_old + ypnew * weight_new;
        const amrex::ParticleReal zp = z_old_p * weight_old + z_new_p * weight_new;
        const amrex::ParticleReal uxp = m_uxpold[i_src] * weight_old
                                      + m_uxpnew[i_src] * weight_new;
        const amrex::ParticleReal uyp = m_uypold[i_src] * weight_old
                                      + m_uypnew[i_src] * weight_new;
        const amrex::ParticleReal uzp = uz_old_p * weight_old
                                      + uz_new_p * weight_new;
        dst.m_aos[i_dst].pos(0) = xp;
#if (AMREX_SPACEDIM == 3)
        dst.m_aos[i_dst].pos(1) = yp;
        dst.m_aos[i_dst].pos(2) = zp;
#elif (AMREX_SPACEDIM == 2)
        dst.m_aos[i_dst].pos(1) = zp;
        amrex::ignore_unused(yp);
#endif
        dst.m_rdata[PIdx::w][i_dst] = m_wpnew[i_src];
        dst.m_rdata[PIdx::ux][i_dst] = uxp;
        dst.m_rdata[PIdx::uy][i_dst] = uyp;
        dst.m_rdata[PIdx::uz][i_dst] = uzp;
    }

    GetParticlePosition m_get_position;

    amrex::ParticleReal* AMREX_RESTRICT m_xpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_ypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_zpold = nullptr;

    amrex::ParticleReal* AMREX_RESTRICT m_uxpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_uypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT m_uzpold = nullptr;

    const amrex::ParticleReal* AMREX_RESTRICT m_uxpnew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uypnew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uzpnew = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_wpnew = nullptr;

    amrex::Real m_gammaboost;
    amrex::Real m_betaboost;
    amrex::Real m_Phys_c;
    amrex::Real m_uzfrm;
    amrex::Real m_inv_c2;
    amrex::Real m_t_boost;
    amrex::Real m_dt;
    amrex::Real m_t_lab;
};

/**
 * \brief Functor to compute Lorentz Transform and store the selected particles in existing
 * particle buffers
 */
class
BackTransformParticleFunctor final : public ComputeParticleDiagFunctor
{
public:
    BackTransformParticleFunctor(WarpXParticleContainer *pc_src, std::string species_name, int num_buffers);
    //void operator ()() const override;
    void operator () (ParticleContainer& pc_dst, int &TotalParticleCounter, int i_buffer) const override;
    void InitData() override;
    void PrepareFunctorData ( int i_buffer, bool ZSliceInDomain, amrex::Real old_z_boost,
                              amrex::Real current_z_boost, amrex::Real t_lab) override;
private:
    WarpXParticleContainer* m_pc_src = nullptr;
    std::string m_species_name;
    int m_num_buffers;
    amrex::Vector<amrex::Real> m_current_z_boost;
    amrex::Vector<amrex::Real> m_old_z_boost;
    amrex::Vector<amrex::Real> m_t_lab;
    amrex::Vector<int> m_perform_backtransform;
};




#endif //
