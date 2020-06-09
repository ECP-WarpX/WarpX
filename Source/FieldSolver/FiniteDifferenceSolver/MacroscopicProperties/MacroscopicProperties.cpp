#include "MacroscopicProperties.H"
#include <AMReX_ParmParse.H>

using namespace amrex;

MacroscopicProperties::MacroscopicProperties ()
{
    ReadParameters();
}

void
MacroscopicProperties::ReadParameters ()
{
    ParmParse pp("macroscopic");
    // Since macroscopic maxwell solve is turned on, user must define sigma, mu, and epsilon //
    pp.get("sigma", m_sigma);
    pp.get("mu", m_mu);
    pp.get("epsilon", m_epsilon);

}

