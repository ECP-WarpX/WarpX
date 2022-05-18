/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "VelocityProperties.H"

#include "Utils/TextMsg.H"

VelocityProperties::VelocityProperties (amrex::ParmParse& pp) {
    // Set defaults
    std::string vel_dist_s = "constant";
    std::string vel_dir_s = "x";
    m_velocity = 0;

    pp.query("bulk_vel_dir", vel_dir_s);
    if(vel_dir_s[0] == '-'){
        m_sign_dir = -1;
    }
    else {
        m_sign_dir = 1;
    }

    if ((vel_dir_s == "x" || vel_dir_s[1] == 'x') ||
       (vel_dir_s == "X" || vel_dir_s[1] == 'X')){
        m_dir = 0;
    }
    else if ((vel_dir_s == "y" || vel_dir_s[1] == 'y') ||
               (vel_dir_s == "Y" || vel_dir_s[1] == 'Y')){
        m_dir = 1;
    }
    else if ((vel_dir_s == "z" || vel_dir_s[1] == 'z') ||
            (vel_dir_s == "Z" || vel_dir_s[1] == 'Z')) {
        m_dir = 2;
    }
    else {
        amrex::Abort(Utils::TextMsg::Err(
            "Cannot interpret <s_name>.bulk_vel_dir input '" + vel_dir_s +
            "'. Please enter +/- x, y, or z with no whitespace between the sign and"+
            " other character."));
    }

    pp.query("beta_distribution_type", vel_dist_s);
    if (vel_dist_s == "constant") {
        queryWithParser(pp, "beta", m_velocity);
        m_type = VelConstantValue;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_velocity > -1 && m_velocity < 1,
            "Magnitude of velocity beta = " + std::to_string(m_velocity) +
            " is greater than or equal to 1"
        );
    }
    else if (vel_dist_s == "parser") {
        std::string str_beta_function;
        Store_parserString(pp, "beta_function(x,y,z)", str_beta_function);
        m_ptr_velocity_parser =
            std::make_unique<amrex::Parser>(makeParser(str_beta_function,{"x","y","z"}));
        m_type = VelParserFunction;
    }
    else {
        amrex::Abort(Utils::TextMsg::Err(
            "Velocity distribution type '" + vel_dist_s + "' not recognized."));
    }
}
