#include "ParticleDiag.H"
#include <AMReX_ParmParse.H>

using namespace amrex;

ParticleDiag::ParticleDiag(std::string diag_name, std::string name, WarpXParticleContainer* pc)
    : m_diag_name(diag_name), m_name(name), m_pc(pc)
{
    ParmParse pp(diag_name + "." + name);
    if (!pp.queryarr("variables", variables)){
        variables = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "ux", "uy", "uz", "w"};
    }
}
