/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "TemperatureProperties.H"

/*
 * Construct TemperatureProperties based on the passed parameters.
 * If temperature is a constant, store value. If a parser, make and
 * store the parser function
 */
TemperatureProperties::TemperatureProperties (amrex::ParmParse& pp) {
    // Set defaults
    amrex::Real theta;
    std::string temp_dist_s = "constant";
    std::string mom_dist_s;

    pp.query("theta_distribution_type", temp_dist_s);
    pp.query("momentum_distribution_type", mom_dist_s);
    if (temp_dist_s == "constant") {
        if (!queryWithParser(pp, "theta", theta)) {
            std::string err_str =  "Temperature parameter theta not specified";
            amrex::Abort(err_str);
        }
        // Do validation on theta value
        if (theta < 0) {
            std::stringstream stringstream;
            stringstream << "Temperature parameter theta = " << theta <<
                " is less than zero, which is not allowed";
            amrex::Abort(stringstream.str().c_str());
        }
        if (mom_dist_s == "maxwell_boltzmann" && theta > 0.01) {
            std::stringstream warnstream;
            warnstream << " Warning: Maxwell-Boltzmann distribution has errors greater than 1%"
                << " for temperature parameter theta > 0.01. (theta = " << theta << " given).";
            amrex::Warning(warnstream.str());
        }
        else if (mom_dist_s == "maxwell_juttner" && theta < 0.1) {
            std::stringstream stringstream;
            stringstream << "Temperature parameter theta = " << theta <<
                " is less than minimum 0.1 allowed for Maxwell-Juttner.";
            amrex::Abort(stringstream.str().c_str());
        }

        m_type = TempConstantValue;
        m_temperature = theta;
    }
    else if (temp_dist_s == "parser") {
        std::string str_theta_function;
        Store_parserString(pp, "theta_function(x,y,z)", str_theta_function);
        m_ptr_temperature_parser =
            std::make_unique<amrex::Parser>(makeParser(str_theta_function,{"x","y","z"}));
        m_type = TempParserFunction;
    }
    else {
        std::stringstream stringstream;
        std::string string;
        stringstream << "Temperature distribution type '" << temp_dist_s << "' not recognized.";
        string = stringstream.str();
        amrex::Abort(string.c_str());
    }
}
