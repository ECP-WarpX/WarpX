#ifndef WARPX_BACKTRANSFORMPARTICLEFUNCTOR_H_
#define WARPX_BACKTRANSFORMPARTICLEFUNCTOR_H_

#include "ComputeParticleDiagFunctor.H"

#include <AMReX.H>
#include <AMReX_AmrParticles.H>

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
    void operator () (ParticleContainer& pc_dst, int i_buffer) const override;
    void InitData() override;
    void PrepareFunctorData ( int i_buffer, bool ZSliceInDomain, amrex::Real old_z_boost,
                              amrex::Real current_z_boost, amrex::Real t_lab) override;
private:
    int m_num_buffers;
    std::string m_species_name;
    WarpXParticleContainer* m_pc_src = nullptr;
    amrex::Vector<amrex::Real> m_current_z_boost;
    amrex::Vector<amrex::Real> m_old_z_boost;
    amrex::Vector<amrex::Real> m_t_lab;
    amrex::Vector<int> m_perform_backtransform;
};

#endif //
