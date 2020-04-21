#ifndef WARPX_PARTICLEDIAG_H_
#define WARPX_PARTICLEDIAG_H_

#include "Particles/WarpXParticleContainer.H"

class ParticleDiag
{
public:
    ParticleDiag(std::string diag_name, std::string name, WarpXParticleContainer* pc);
private:
    std::string m_diag_name;
    std::string m_name;
    amrex::Vector< std::string > variables;
    WarpXParticleContainer* m_pc;
};
    
#endif // WARPX_PARTICLEDIAG_H_
