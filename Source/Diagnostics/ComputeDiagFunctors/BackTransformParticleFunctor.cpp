#include "BackTransformParticleFunctor.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "WarpX.H"
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BaseFwd.H>

/**
 * \brief Functor to compute Lorentz Transform and store the selected particles in existing
 * particle buffers
 */
BackTransformParticleFunctor::BackTransformParticleFunctor (
                              WarpXParticleContainer *pc_src,
                              std::string species_name,
                              int num_buffers)
    : m_pc_src(pc_src), m_species_name(species_name), m_num_buffers(num_buffers)
{
    amrex::Print() << " in functor : sp name: " << species_name << "\n";;
    InitData();
}

struct SelectParticles
{
    using TmpParticles = WarpXParticleContainer::TmpParticles;

    SelectParticles( const WarpXParIter& a_pti, TmpParticles& tmp_particle_data,
                     amrex::Real current_z_boost, amrex::Real old_z_boost,
                     int a_offset = 0)
        : m_current_z_boost(current_z_boost), m_old_z_boost(old_z_boost)
    {
        if (tmp_particle_data.size() == 0) return;
        m_get_position = GetParticlePosition(a_pti, a_offset);
 
        const auto lev = a_pti.GetLevel();
        const auto index = a_pti.GetPairIndex();

        zpold = tmp_particle_data[lev][index][TmpIdx::zold].dataPtr();
    }

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
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
                                amrex::Real t_lab, int a_offset = 0)
        : m_t_boost(t_boost), m_dt(dt), m_t_lab(t_lab)
    {
        using namespace amrex::literals;
        
        if (tmp_particle_data.size() == 0) return;
        m_get_position = GetParticlePosition(a_pti, a_offset);
        
        auto& attribs = a_pti.GetAttribs();
        m_wpnew = attribs[PIdx::w].dataPtr(); 
        m_uxpnew = attribs[PIdx::ux].dataPtr(); 
        m_uypnew = attribs[PIdx::ux].dataPtr(); 
        m_uzpnew = attribs[PIdx::ux].dataPtr(); 

        const auto lev = a_pti.GetLevel();
        const auto index = a_pti.GetPairIndex();

        m_xpold = tmp_particle_data[lev][index][TmpIdx::xold].dataPtr();
        m_ypold = tmp_particle_data[lev][index][TmpIdx::yold].dataPtr();
        m_zpold = tmp_particle_data[lev][index][TmpIdx::zold].dataPtr();
        m_uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
        m_uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
        m_uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();

        m_betaboost = WarpX::beta_boost;
        m_gammaboost = WarpX::gamma_boost;
        m_Phys_c = PhysConst::c;
        m_inv_c2 = 1._rt/(m_Phys_c * m_Phys_c);
        m_uzfrm = -m_gammaboost*m_betaboost*m_Phys_c;
    }

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator () (const DstData& dst, const SrcData& src, int i_src, int i_dst) const noexcept
    {
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
        const amrex::ParticleReal uzp = m_uzpold[i_src] * weight_old
                                      + m_uzpnew[i_src] * weight_new;
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

void
BackTransformParticleFunctor::operator () (ParticleContainer& pc_dst, int i_buffer) const
{
    amrex::Print() << " in BTD functor operator \n";
    ParticleContainer pc_tmp(&WarpX::GetInstance());
    auto &warpx = WarpX::GetInstance();
    // get particle slice
    amrex::Print() << " finest level : " << m_pc_src->finestLevel() << "\n"; 
    const int nlevs = std::max(0, m_pc_src->finestLevel()+1);
    amrex::Print() << " nlevs : " << nlevs << "\n";
    auto tmp_particle_data = m_pc_src->getTmpParticleData();

    for (int lev = 0; lev < nlevs; ++lev) {
        amrex::Real t_boost = warpx.gett_new(0);
        amrex::Real dt = warpx.getdt(0);

        for (WarpXParIter pti(*m_pc_src, lev); pti.isValid(); ++pti) {
            auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
            auto ptile_tmp = pc_tmp.DefineAndReturnParticleTile(lev, pti.index(), pti.LocalTileIndex() );
        }

        auto& particles = m_pc_src->GetParticles(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            // Temporary arrays to store copy_flag and copy_index for particles
            // that cross the z-slice
            amrex::Gpu::DeviceVector<int> FlagForPartCopy;
            amrex::Gpu::DeviceVector<int> IndexForPartCopy;
            
            for (WarpXParIter pti(*m_pc_src, lev); pti.isValid(); ++pti) {

                auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

                const auto GetPosition = GetParticlePosition(pti);
                const auto GetParticleFilter = SelectParticles(pti, tmp_particle_data,
                                               m_current_z_boost[i_buffer],
                                               m_old_z_boost[i_buffer]);
                const auto GetParticleLorentzTransform = LorentzTransformParticles(
                                                         pti, tmp_particle_data,
                                                         t_boost, dt,
                                                         m_t_lab[i_buffer]);

                long const np = pti.numParticles();

                FlagForPartCopy.resize(np);
                IndexForPartCopy.resize(np);

                int* const AMREX_RESTRICT Flag = FlagForPartCopy.dataPtr();
                int* const AMREX_RESTRICT IndexLocation = IndexForPartCopy.dataPtr();
                
                auto& ptile_src = particles[index];
                // Flag particles that need to be copied if they cross the z-slice
                // setting this to 1 for testing (temporarily)
                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE(int i)
                {
                    Flag[i] = GetParticleFilter(ptile_src, i);
                });

                const int total_partdiag_size = amrex::Scan::ExclusiveSum(np,Flag,IndexLocation);
                auto& ptile_tmp = pc_tmp.DefineAndReturnParticleTile(lev, pti.index(), pti.LocalTileIndex() );
                auto old_size = ptile_tmp.numParticles();
                ptile_tmp.resize(old_size + total_partdiag_size); 
                amrex::Print() << " total partdiag_size : " << total_partdiag_size << "\n";
                amrex::Print() << " old size : " << old_size << "\n";
                amrex::Print() << " total part diag size :  " << total_partdiag_size << "\n";

                auto count = amrex::filterParticles(ptile_tmp, ptile_src, GetParticleFilter, 0, old_size, np);
                amrex::Print() << " count : " << count << "\n";
                auto dst_data = ptile_tmp.getParticleTileData();
                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE(int i)
                {
                   if (Flag[i] == 1) GetParticleLorentzTransform(dst_data, ptile_src, i,
                                                                 old_size + IndexLocation[i]);
                });
                //auto transformedParticles = amrex::filterAndTransformParticles(ptile_tmp,
                //                            ptile_src, FlagForPartCopy.dataPtr(),
                //                            GetParticleLorentzTransform, 0, old_size, np);

            }
        }        
        amrex::Gpu::synchronize();        
    }
}


void
BackTransformParticleFunctor::InitData()
{
    m_current_z_boost.resize(m_num_buffers);
    m_old_z_boost.resize(m_num_buffers);
    m_t_lab.resize(m_num_buffers);
}

void
BackTransformParticleFunctor::PrepareFunctorData ( int i_buffer, amrex::Real old_z_boost,
                              amrex::Real current_z_boost, amrex::Real t_lab)
{
    m_old_z_boost[i_buffer] = old_z_boost;
    m_current_z_boost[i_buffer] = current_z_boost;
    m_t_lab[i_buffer] = t_lab;
}
