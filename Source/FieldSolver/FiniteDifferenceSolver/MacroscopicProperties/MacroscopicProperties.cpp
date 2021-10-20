#include "MacroscopicProperties.H"

#include "Utils/WarpXUtil.H"
#include "WarpX.H"

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

#include <AMReX_BaseFwd.H>

#include <memory>
#include <sstream>

using namespace amrex;

GetSigmaMacroparameter::GetSigmaMacroparameter () noexcept
{
    auto& warpx = WarpX::GetInstance();
    auto& macroscopic_properties = warpx.GetMacroscopicProperties();
    if (macroscopic_properties.m_sigma_s == "constant") {
        m_type = ConstantValue;
        m_value = macroscopic_properties.m_sigma;
    }
    else if (macroscopic_properties.m_sigma_s == "parse_sigma_function") {
        m_type = ParserFunction;
        m_parser = macroscopic_properties.m_sigma_parser->compile<3>();
    }
}

GetMuMacroparameter::GetMuMacroparameter () noexcept
{
    auto& warpx = WarpX::GetInstance();
    auto& macroscopic_properties = warpx.GetMacroscopicProperties();
    if (macroscopic_properties.m_mu_s == "constant") {
        m_type = ConstantValue;
        m_value = macroscopic_properties.m_mu;
    }
    else if (macroscopic_properties.m_mu_s == "parse_mu_function") {
        m_type = ParserFunction;
        m_parser = macroscopic_properties.m_mu_parser->compile<3>();
    }
}

GetEpsilonMacroparameter::GetEpsilonMacroparameter () noexcept
{
    auto& warpx = WarpX::GetInstance();
    auto& macroscopic_properties = warpx.GetMacroscopicProperties();
    if (macroscopic_properties.m_epsilon_s == "constant") {
        m_type = ConstantValue;
        m_value = macroscopic_properties.m_epsilon;
    }
    else if (macroscopic_properties.m_epsilon_s == "parse_epsilon_function") {
        m_type = ParserFunction;
        m_parser = macroscopic_properties.m_epsilon_parser->compile<3>();
    }
}


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
    if (queryWithParser(pp_macroscopic, "sigma", m_sigma)) {
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
        WarpX::GetInstance().RecordWarning("Macroscopic properties",
            warnMsg.str());
    }
    // initialization of sigma (conductivity) with parser
    if (m_sigma_s == "parse_sigma_function") {
        Store_parserString(pp_macroscopic, "sigma_function(x,y,z)", m_str_sigma_function);
        m_sigma_parser = std::make_unique<Parser>(
                                 makeParser(m_str_sigma_function,{"x","y","z"}));
    }

    bool epsilon_specified = false;
    if (queryWithParser(pp_macroscopic, "epsilon", m_epsilon)) {
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
        WarpX::GetInstance().RecordWarning("Macroscopic properties",
            warnMsg.str());
    }

    // initialization of epsilon (permittivity) with parser
    if (m_epsilon_s == "parse_epsilon_function") {
        Store_parserString(pp_macroscopic, "epsilon_function(x,y,z)", m_str_epsilon_function);
        m_epsilon_parser = std::make_unique<Parser>(
                                 makeParser(m_str_epsilon_function,{"x","y","z"}));
    }

    // Query input for material permittivity, epsilon.
    bool mu_specified = false;
    if (queryWithParser(pp_macroscopic, "mu", m_mu)) {
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
        WarpX::GetInstance().RecordWarning("Macroscopic properties",
            warnMsg.str());
    }

    // initialization of mu (permeability) with parser
    if (m_mu_s == "parse_mu_function") {
        Store_parserString(pp_macroscopic, "mu_function(x,y,z)", m_str_mu_function);
        m_mu_parser = std::make_unique<Parser>(
                                 makeParser(m_str_mu_function,{"x","y","z"}));
    }

}

void
MacroscopicProperties::InitData ()
{
    amrex::Print() << "we are in init data of macro \n";
    auto & warpx = WarpX::GetInstance();

    IntVect Ex_stag = warpx.getEfield_fp(0,0).ixType().toIntVect();
    IntVect Ey_stag = warpx.getEfield_fp(0,1).ixType().toIntVect();
    IntVect Ez_stag = warpx.getEfield_fp(0,2).ixType().toIntVect();
    IntVect Bx_stag = warpx.getBfield_fp(0,0).ixType().toIntVect();
    IntVect By_stag = warpx.getBfield_fp(0,1).ixType().toIntVect();
    IntVect Bz_stag = warpx.getBfield_fp(0,2).ixType().toIntVect();


    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Ex_IndexType[idim]      = Ex_stag[idim];
        Ey_IndexType[idim]      = Ey_stag[idim];
        Ez_IndexType[idim]      = Ez_stag[idim];
        Bx_IndexType[idim]      = Bx_stag[idim];
        By_IndexType[idim]      = By_stag[idim];
        Bz_IndexType[idim]      = Bz_stag[idim];
    }
#if (AMREX_SPACEDIM==2)
        Ex_IndexType[2]      = 0;
        Ey_IndexType[2]      = 0;
        Ez_IndexType[2]      = 0;
        Bx_IndexType[2]      = 0;
        By_IndexType[2]      = 0;
        Bz_IndexType[2]      = 0;
#endif


}

