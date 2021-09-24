/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <TemperatureProperties.H>

/* 
 * Construct TemperatureProperties based on the passed parameters.
 * If temperature is a constant, store value. If a parser, make and
 * store the parser function
 */
TemperatureProperties::TemperatureProperties (amrex::ParmParse& pp) {
    // Set defaults
    amrex::Real theta = 20.;
    std::string temp_dist_s = "constant";

    pp.query("theta_distribution_type", temp_dist_s);
    if (temp_dist_s == "constant") {
        queryWithParser(pp, "theta", theta);
        m_type = ConstantValue;
        m_temperature = theta;
    }
    else if (temp_dist_s == "parser") {
        std::string str_theta_function;
        Store_parserString(pp, "theta_function(x,y,z)", str_theta_function);
        m_ptr_temperature_parser =
            std::make_unique<amrex::Parser>(makeParser(str_theta_function,{"x","y","z"}));
        m_type = ParserFunction;
    }
    else {
        std::stringstream stringstream;
        std::string string;
        stringstream << "Temperature distribution type '" << temp_dist_s << "' not recognized." << std::endl;
        string = stringstream.str();
        amrex::Abort(string.c_str());
    }
}
