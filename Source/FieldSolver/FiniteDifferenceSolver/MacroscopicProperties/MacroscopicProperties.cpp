#include "MacroscopicProperties.H"
#include <AMReX_ParmParse.H>
#include "WarpX.H"

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
//    pp.get("sigma", m_sigma);
//    pp.get("mu", m_mu);
//    pp.get("epsilon", m_epsilon);

    pp.get("sigma_init_style", m_sigma_s);
    if (m_sigma_s == "constant") pp.get("sigma", m_sigma);

    pp.get("epsilon_init_style", m_epsilon_s);
    if (m_epsilon_s == "constant") pp.get("epsilon", m_epsilon);

    pp.get("mu_init_style", m_mu_s);
    if (m_mu_s == "constant") pp.get("mu", m_mu);

}

void
MacroscopicProperties::InitData ()
{
    amrex::Print() << "we are in init data of macro \n";
    auto & warpx = WarpX::GetInstance();

    // Get BoxArray and DistributionMap of warpX.
    int lev = 0;
    BoxArray ba = warpx.boxArray(lev);
    DistributionMapping dmap = warpx.DistributionMap(lev);
    int ng = 1;
      // sigma
    m_sigma_mf = std::make_unique<MultiFab>(ba, dmap, 1, ng); // cell-centered
    // for now initializing with constant sigma
    m_sigma_mf->setVal(m_sigma);
      // eps - cell-centered
    m_eps_mf = std::make_unique<MultiFab>(ba, dmap, 1, ng);
    m_eps_mf->setVal(m_epsilon);
      // mu - node-based
    m_mu_mf = std::make_unique<MultiFab>(amrex::convert(ba,amrex::IntVect::TheUnitVector()), dmap, 1, ng);
    m_mu_mf->setVal(m_mu);

}



