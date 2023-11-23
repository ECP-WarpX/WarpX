/* Copyright 2023 Luca Fedeli
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ExternalField.H"

#include "Utils/TextMsg.H"
#include "Utils/Parser/ParserUtils.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <algorithm>
#include <vector>

namespace
{
   ExternalFieldType string_to_external_field_type(std::string s)
   {
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);

        if ( s == "" || s == "default")
            return ExternalFieldType::default_zero;
        else if ( s == "constant")
            return ExternalFieldType::constant;
        else if ( s == "parse_b_ext_grid_function" || s == "parse_e_ext_grid_function")
            return ExternalFieldType::parse_ext_grid_function;
        else if ( s == "read_from_file")
            return ExternalFieldType::read_from_file;
        else
            WARPX_ABORT_WITH_MESSAGE(
                "'" + s + "' is an unknown external field type!");

        return ExternalFieldType::default_zero;
   }
}

std::unique_ptr<ExternalFieldParams> ReadExternalFieldParams(const amrex::ParmParse& pp_warpx)
{
    auto p_external_field = std::make_unique<ExternalFieldParams>();

    // default values of E_external_grid and B_external_grid
    // are used to set the E and B field when "constant" or
    // "parser" is not explicitly used in the input.
    std::string B_ext_grid_s = "";
    pp_warpx.query("B_ext_grid_init_style", B_ext_grid_s);
    p_external_field->B_ext_grid_type = string_to_external_field_type(B_ext_grid_s);

    std::string E_ext_grid_s = "";
    pp_warpx.query("E_ext_grid_init_style", E_ext_grid_s);
    p_external_field->E_ext_grid_type = string_to_external_field_type(E_ext_grid_s);

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    auto v_B = std::vector<amrex::Real>(3);
    if (p_external_field->B_ext_grid_type == ExternalFieldType::constant)
        utils::parser::getArrWithParser(pp_warpx, "B_external_grid", v_B);
    std::copy(v_B.begin(), v_B.end(), p_external_field->B_external_grid.begin());

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    auto v_E = std::vector<amrex::Real>(3);
    if (p_external_field->E_ext_grid_type == ExternalFieldType::constant)
        utils::parser::getArrWithParser(pp_warpx, "E_external_grid", v_E);
    std::copy(v_E.begin(), v_E.end(), p_external_field->E_external_grid.begin());




    // if the input string for the B-field is "parse_b_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (p_external_field->B_ext_grid_type == ExternalFieldType::parse_ext_grid_function) {

        //! Strings storing parser function to initialize the components of the magnetic field on the grid
        std::string str_Bx_ext_grid_function;
        std::string str_By_ext_grid_function;
        std::string str_Bz_ext_grid_function;

#ifdef WARPX_DIM_RZ
        std::stringstream warnMsg;
        warnMsg << "Parser for external B (r and theta) fields does not work with RZ\n"
            << "The initial Br and Bt fields are currently hardcoded to 0.\n"
            << "The initial Bz field should only be a function of z.\n";
        ablastr::warn_manager::WMRecordWarning(
          "Inputs", warnMsg.str(), ablastr::warn_manager::WarnPriority::high);
        str_Bx_ext_grid_function = "0";
        str_By_ext_grid_function = "0";
#else
        utils::parser::Store_parserString(pp_warpx, "Bx_external_grid_function(x,y,z)",
          str_Bx_ext_grid_function);
        utils::parser::Store_parserString(pp_warpx, "By_external_grid_function(x,y,z)",
          str_By_ext_grid_function);
#endif
        utils::parser::Store_parserString(pp_warpx, "Bz_external_grid_function(x,y,z)",
            str_Bz_ext_grid_function);

        p_external_field->Bxfield_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(str_Bx_ext_grid_function,{"x","y","z"}));
        p_external_field->Byfield_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(str_By_ext_grid_function,{"x","y","z"}));
        p_external_field->Bzfield_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(str_Bz_ext_grid_function,{"x","y","z"}));
    }



    // if the input string for the E-field is "parse_e_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (p_external_field->E_ext_grid_type == ExternalFieldType::parse_ext_grid_function) {

#ifdef WARPX_DIM_RZ
        WARPX_ABORT_WITH_MESSAGE(
            "E parser for external fields does not work with RZ -- TO DO");
#endif

        //! Strings storing parser function to initialize the components of the electric field on the grid
        std::string str_Ex_ext_grid_function;
        std::string str_Ey_ext_grid_function;
        std::string str_Ez_ext_grid_function;

        utils::parser::Store_parserString(pp_warpx, "Ex_external_grid_function(x,y,z)",
            str_Ex_ext_grid_function);
        utils::parser::Store_parserString(pp_warpx, "Ey_external_grid_function(x,y,z)",
           str_Ey_ext_grid_function);
        utils::parser::Store_parserString(pp_warpx, "Ez_external_grid_function(x,y,z)",
           str_Ez_ext_grid_function);

        p_external_field->Exfield_parser = std::make_unique<amrex::Parser>(
           utils::parser::makeParser(str_Ex_ext_grid_function,{"x","y","z"}));
        p_external_field->Eyfield_parser = std::make_unique<amrex::Parser>(
           utils::parser::makeParser(str_Ey_ext_grid_function,{"x","y","z"}));
        p_external_field->Ezfield_parser = std::make_unique<amrex::Parser>(
           utils::parser::makeParser(str_Ez_ext_grid_function,{"x","y","z"}));
    }


    if (p_external_field->E_ext_grid_type == ExternalFieldType::read_from_file ||
        p_external_field->B_ext_grid_type == ExternalFieldType::read_from_file){
            std::string read_fields_from_path="./";
            pp_warpx.query("read_fields_from_path", p_external_field->external_fields_path);
    }

    return p_external_field;
}
