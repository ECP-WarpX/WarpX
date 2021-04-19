#include "BackTransformParticleFunctor.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "WarpX.H"
#include <AMReX.H>
#include <AMReX_Print.H>

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

        xpold = tmp_particle_data[lev][index][TmpIdx::xold].dataPtr();
        ypold = tmp_particle_data[lev][index][TmpIdx::yold].dataPtr();
        zpold = tmp_particle_data[lev][index][TmpIdx::zold].dataPtr();
        uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
        uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
        uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();
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
    
    amrex::ParticleReal* AMREX_RESTRICT xpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT ypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT zpold = nullptr;

    amrex::ParticleReal* AMREX_RESTRICT uxpold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT uypold = nullptr;
    amrex::ParticleReal* AMREX_RESTRICT uzpold = nullptr;
};


void
BackTransformParticleFunctor::operator () (ParticleContainer& pc_dst, int i_buffer) const
{
    amrex::Print() << " in BTD functor operator \n";
    ParticleContainer pc_tmp(&WarpX::GetInstance());
    // get particle slice
    amrex::Print() << " finest level : " << m_pc_src->finestLevel() << "\n"; 
    const int nlevs = std::max(0, m_pc_src->finestLevel()+1);
    amrex::Print() << " nlevs : " << nlevs << "\n";
    auto tmp_particle_data = m_pc_src->getTmpParticleData();
    for (int lev = 0; lev < nlevs; ++lev) {
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
                auto ptile_tmp = pc_tmp.DefineAndReturnParticleTile(lev, pti.index(), pti.LocalTileIndex() );
                auto old_size = ptile_tmp.numParticles();
                ptile_tmp.resize(old_size + total_partdiag_size); 
                amrex::Print() << " total partdiag_size : " << total_partdiag_size << "\n";
                amrex::Print() << " old size : " << old_size << "\n";
                amrex::Print() << " total part diag size :  " << total_partdiag_size << "\n";

                auto count = amrex::filterParticles(ptile_tmp, ptile_src, GetParticleFilter, 0, old_size, np);
                amrex::Print() << " count : " << count << "\n";

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
