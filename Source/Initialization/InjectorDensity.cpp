/* Copyright 2019-2020 Axel Huebl, Ligia Diana Amorim, Maxence Thevenet
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "InjectorDensity.H"
#include "WarpX.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <cctype>
#include <vector>

using namespace amrex;

void InjectorDensity::clear ()
{
    switch (type)
    {
        case Type::parser:
        {
            break;
        }
        case Type::predefined:
        {
            object.predefined.clear();
            break;
        }
        case Type::fromfile:
        {
            object.fromfile.clear();
            break;
        }
        default:
            return;
    }
}

InjectorDensityPredefined::InjectorDensityPredefined (
        std::string const& a_species_name) noexcept
{
    const ParmParse pp_species_name(a_species_name);

    std::vector<amrex::Real> v;
    // Read parameters for the predefined plasma profile.
    utils::parser::getArrWithParser(
            pp_species_name, "predefined_profile_params", v);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v.size() <= 6,
                                     "Too many parameters for InjectorDensityPredefined");
    for (int i = 0; i < static_cast<int>(v.size()); ++i) {
        p[i] = v[i];
    }

    // Parse predefined profile name, and update member variable profile.
    std::string which_profile_s;
    pp_species_name.query("predefined_profile_name", which_profile_s);
    std::transform(which_profile_s.begin(), which_profile_s.end(),
                   which_profile_s.begin(), ::tolower);
    if (which_profile_s == "parabolic_channel"){
        profile = Profile::parabolic_channel;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v.size() > 5,
                                         "InjectorDensityPredefined::parabolic_channel: not enough parameters");
    }
}

InjectorDensityFromFile::InjectorDensityFromFile (std::string const & a_species_name)
{
    //get the file path for reading the external density and the name of the field to be read
    const ParmParse pp (a_species_name);
    const ParmParse pp_warpx ("warpx");
    const ParmParse pp_amr ("amr");

    std::string external_density_path, density;
    pp.query("read_density_from_path", external_density_path);
    pp.query("density_name", density);

    // Get WarpX domain info
    WarpX& warpx = WarpX::GetInstance();
    amrex::Geometry const& geom0 = warpx.Geom(0);
    m_cell_size = geom0.CellSizeArray();
    auto geom_lo = geom0.ProbLoArray();
    auto geom_hi = geom0.ProbHiArray();
    lo0 = geom_lo[0];
    lo1 = geom_lo[1];
    lo2 = geom_lo[2];
    hi0 = geom_hi[0];
    hi1 = geom_hi[1];
    hi2 = geom_hi[2];

    // creating the mulitfab array
    m_rho = warpx.getFieldPointer(warpx::fields::FieldType::rho_fp, 0);
    m_rho->setVal(0);

    WarpX::ReadExternalFieldFromFile( external_density_path, m_rho, density, "" );

}

// Note that we are not allowed to have non-trivial destructor.
// So we rely on clear() to free memory if needed.
void InjectorDensityPredefined::clear ()
{
}
void InjectorDensityFromFile::clear ()
{
}
