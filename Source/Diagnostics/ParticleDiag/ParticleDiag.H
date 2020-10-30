#ifndef WARPX_PARTICLEDIAG_H_
#define WARPX_PARTICLEDIAG_H_

#include "Parser/WarpXParserWrapper.H"
#include "Utils/WarpXUtil.H"
#include "Particles/WarpXParticleContainer.H"

#include <memory>

class ParticleDiag
{
public:
    ParticleDiag(std::string diag_name, std::string name, WarpXParticleContainer* pc);
    WarpXParticleContainer* getParticleContainer() const { return m_pc; }
    std::string getSpeciesName() const { return m_name; }
    amrex::Vector<int> plot_flags;

    bool m_do_random_filter  = false;
    bool m_do_uniform_filter = false;
    bool m_do_parser_filter  = false;
    bool m_do_geom_filter    = false;
    amrex::Real m_random_fraction = 1.0;
    int m_uniform_stride = 1;
    std::unique_ptr<ParserWrapper<7>> m_particle_filter_parser;
    amrex::RealBox m_diag_domain;

private:
    std::string m_diag_name;
    std::string m_name;
    amrex::Vector< std::string > variables;
    WarpXParticleContainer* m_pc;
};

#endif // WARPX_PARTICLEDIAG_H_
