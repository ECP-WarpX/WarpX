#include "MacroscopicProperties.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_Array4.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Parser.H>

#include <AMReX_BaseFwd.H>

#include <memory>
#include <sstream>

using namespace amrex;

MacroscopicProperties::MacroscopicProperties ()
{
    ReadParameters();
}

void
MacroscopicProperties::ReadParameters ()
{
    ParmParse pp_macroscopic("macroscopic");
    // Since macroscopic maxwell solve is turned on,
    // user-defined sigma, mu, and epsilon are queried.
    // The vacuum values are used as default for the macroscopic parameters
    // with a warning message to the user to indicate that no value was specified.

    // Query input for material conductivity, sigma.
    bool sigma_specified = false;
    if (utils::parser::queryWithParser(pp_macroscopic, "sigma", m_sigma)) {
        m_sigma_s = "constant";
        sigma_specified = true;
    }
    if (pp_macroscopic.query("sigma_function(x,y,z)", m_str_sigma_function) ) {
        m_sigma_s = "parse_sigma_function";
        sigma_specified = true;
    }
    if (!sigma_specified) {
        std::stringstream warnMsg;
        warnMsg << "Material conductivity is not specified. Using default vacuum value of " <<
            m_sigma << " in the simulation.";
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            warnMsg.str());
    }
    // initialization of sigma (conductivity) with parser
    if (m_sigma_s == "parse_sigma_function") {
        utils::parser::Store_parserString(
            pp_macroscopic, "sigma_function(x,y,z)", m_str_sigma_function);
        m_sigma_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(m_str_sigma_function,{"x","y","z"}));
    }

    bool epsilon_specified = false;
    if (utils::parser::queryWithParser(pp_macroscopic, "epsilon", m_epsilon)) {
        m_epsilon_s = "constant";
        epsilon_specified = true;
    }
    if (pp_macroscopic.query("epsilon_function(x,y,z)", m_str_epsilon_function) ) {
        m_epsilon_s = "parse_epsilon_function";
        epsilon_specified = true;
    }
    if (!epsilon_specified) {
        std::stringstream warnMsg;
        warnMsg << "Material permittivity is not specified. Using default vacuum value of " <<
            m_epsilon << " in the simulation.";
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            warnMsg.str());
    }

    // initialization of epsilon (permittivity) with parser
    if (m_epsilon_s == "parse_epsilon_function") {
        utils::parser::Store_parserString(
            pp_macroscopic, "epsilon_function(x,y,z)", m_str_epsilon_function);
        m_epsilon_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(m_str_epsilon_function,{"x","y","z"}));
    }

    // Query input for material permittivity, epsilon.
    bool mu_specified = false;
    if (utils::parser::queryWithParser(pp_macroscopic, "mu", m_mu)) {
        m_mu_s = "constant";
        mu_specified = true;
    }
    if (pp_macroscopic.query("mu_function(x,y,z)", m_str_mu_function) ) {
        m_mu_s = "parse_mu_function";
        mu_specified = true;
    }
    if (!mu_specified) {
        std::stringstream warnMsg;
        warnMsg << "Material permittivity is not specified. Using default vacuum value of " <<
            m_mu << " in the simulation.";
        ablastr::warn_manager::WMRecordWarning("Macroscopic properties",
            warnMsg.str());
    }

    // initialization of mu (permeability) with parser
    if (m_mu_s == "parse_mu_function") {
        utils::parser::Store_parserString(
            pp_macroscopic, "mu_function(x,y,z)", m_str_mu_function);
        m_mu_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(m_str_mu_function,{"x","y","z"}));
    }

}

void
MacroscopicProperties::InitData ()
{
    amrex::Print() << Utils::TextMsg::Info("we are in init data of macro");
    auto & warpx = WarpX::GetInstance();

    // Get BoxArray and DistributionMap of warpx instance.
    int lev = 0;
    amrex::BoxArray ba = warpx.boxArray(lev);
    amrex::DistributionMapping dmap = warpx.DistributionMap(lev);
    const amrex::IntVect ng_EB_alloc = warpx.getngEB();
    // Define material property multifabs using ba and dmap from WarpX instance
    // sigma is cell-centered MultiFab
    m_sigma_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 1, ng_EB_alloc);
    // epsilon is cell-centered MultiFab
    m_eps_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 1, ng_EB_alloc);
    // mu is cell-centered MultiFab
    m_mu_mf = std::make_unique<amrex::MultiFab>(ba, dmap, 1, ng_EB_alloc);
    // Initialize sigma
    if (m_sigma_s == "constant") {

        m_sigma_mf->setVal(m_sigma);

    } else if (m_sigma_s == "parse_sigma_function") {

        InitializeMacroMultiFabUsingParser(m_sigma_mf.get(), m_sigma_parser->compile<3>(), lev);
    }
    // Initialize epsilon
    if (m_epsilon_s == "constant") {

        m_eps_mf->setVal(m_epsilon);

    } else if (m_epsilon_s == "parse_epsilon_function") {

        InitializeMacroMultiFabUsingParser(m_eps_mf.get(), m_epsilon_parser->compile<3>(), lev);

    }
    // In the Maxwell solver, `epsilon` is used in the denominator.
    // Therefore, it needs to be strictly positive
    bool const local=true;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( m_eps_mf->min(0,0,local) > 0,
    "WarpX encountered zero or negative values for the relative permittivity `epsilon`. Please check the initialization of `epsilon`.");

    // Initialize mu
    if (m_mu_s == "constant") {

        m_mu_mf->setVal(m_mu);

    } else if (m_mu_s == "parse_mu_function") {

        InitializeMacroMultiFabUsingParser(m_mu_mf.get(), m_mu_parser->compile<3>(), lev);

    }

    amrex::IntVect sigma_stag = m_sigma_mf->ixType().toIntVect();
    amrex::IntVect epsilon_stag = m_eps_mf->ixType().toIntVect();
    amrex::IntVect mu_stag = m_mu_mf->ixType().toIntVect();
    amrex::IntVect Ex_stag = warpx.getEfield_fp(0,0).ixType().toIntVect();
    amrex::IntVect Ey_stag = warpx.getEfield_fp(0,1).ixType().toIntVect();
    amrex::IntVect Ez_stag = warpx.getEfield_fp(0,2).ixType().toIntVect();

    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        sigma_IndexType[idim]   = sigma_stag[idim];
        epsilon_IndexType[idim] = epsilon_stag[idim];
        mu_IndexType[idim]      = mu_stag[idim];
        Ex_IndexType[idim]      = Ex_stag[idim];
        Ey_IndexType[idim]      = Ey_stag[idim];
        Ez_IndexType[idim]      = Ez_stag[idim];
        macro_cr_ratio[idim]    = 1;
    }
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        sigma_IndexType[2]   = 0;
        epsilon_IndexType[2] = 0;
        mu_IndexType[2]      = 0;
        Ex_IndexType[2]      = 0;
        Ey_IndexType[2]      = 0;
        Ez_IndexType[2]      = 0;
        macro_cr_ratio[2]    = 1;
#endif
}

void
MacroscopicProperties::InitializeMacroMultiFabUsingParser (
                       amrex::MultiFab *macro_mf,
                       amrex::ParserExecutor<3> const& macro_parser,
                       const int lev)
{
    WarpX& warpx = WarpX::GetInstance();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_lev = warpx.Geom(lev).CellSizeArray();
    const amrex::RealBox& real_box = warpx.Geom(lev).ProbDomain();
    amrex::IntVect iv = macro_mf->ixType().toIntVect();
    for ( amrex::MFIter mfi(*macro_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // Initialize ghost cells in addition to valid cells

        const amrex::Box& tb = mfi.tilebox( iv, macro_mf->nGrowVect());
        amrex::Array4<amrex::Real> const& macro_fab =  macro_mf->array(mfi);
        amrex::ParallelFor (tb,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Shift x, y, z position based on index type
                amrex::Real fac_x = (1._rt - iv[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i * dx_lev[0] + real_box.lo(0) + fac_x;
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real z = j * dx_lev[1] + real_box.lo(1) + fac_z;
#else
                amrex::Real fac_y = (1._rt - iv[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real y = j * dx_lev[1] + real_box.lo(1) + fac_y;
                amrex::Real fac_z = (1._rt - iv[2]) * dx_lev[2] * 0.5_rt;
                amrex::Real z = k * dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // initialize the macroparameter
                macro_fab(i,j,k) = macro_parser(x,y,z);
        });

    }
}
