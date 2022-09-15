/* Copyright 2019-2022 Andrew Myers, Burlen Loring, Luca Fedeli
 * Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParserUtils.H"
#include "Utils/TextMsg.H"

#include <limits>

void utils::parser::Store_parserString(
    const amrex::ParmParse& pp,
    std::string query_string,
    std::string& stored_string)
{
    std::vector<std::string> f;
    pp.getarr(query_string.c_str(), f);
    stored_string.clear();
    for (auto const& s : f) {
        stored_string += s;
    }
    f.clear();
}

int utils::parser::safeCastToInt(
    const amrex::Real x, const std::string& real_name)
{
    int result = 0;
    bool error_detected = false;
    std::string assert_msg;
    // (2.0*(numeric_limits<int>::max()/2+1)) converts numeric_limits<int>::max()+1 to a real ensuring accuracy to all digits
    // This accepts x = 2**31-1 but rejects 2**31.
    using namespace amrex::literals;
    constexpr amrex::Real max_range = (2.0_rt*static_cast<amrex::Real>(std::numeric_limits<int>::max()/2+1));
    if (x < max_range) {
        if (std::ceil(x) >= std::numeric_limits<int>::min()) {
            result = static_cast<int>(x);
        } else {
            error_detected = true;
            assert_msg = "Negative overflow detected when casting " + real_name + " = " + std::to_string(x) + " to int";
        }
    } else if (x > 0) {
        error_detected = true;
        assert_msg =  "Overflow detected when casting " + real_name + " = " + std::to_string(x) + " to int";
    } else {
        error_detected = true;
        assert_msg =  "NaN detected when casting " + real_name + " to int";
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!error_detected, assert_msg);
    return result;
}

Parser makeParser (std::string const& parse_function, amrex::Vector<std::string> const& varnames)
{
    // Since queryWithParser recursively calls this routine, keep track of symbols
    // in case an infinite recursion is found (a symbol's value depending on itself).
    static std::set<std::string> recursive_symbols;

    Parser parser(parse_function);
    parser.registerVariables(varnames);

    std::set<std::string> symbols = parser.symbols();
    for (auto const& v : varnames) symbols.erase(v.c_str());

    // User can provide inputs under this name, through which expressions
    // can be provided for arbitrary variables. PICMI inputs are aware of
    // this convention and use the same prefix as well. This potentially
    // includes variable names that match physical or mathematical
    // constants, in case the user wishes to enforce a different
    // system of units or some form of quasi-physical behavior in the
    // simulation. Thus, this needs to override any built-in
    // constants.
    ParmParse pp_my_constants("my_constants");

    // Physical / Numerical Constants available to parsed expressions
    static std::map<std::string, amrex::Real> warpx_constants =
      {
       {"clight", PhysConst::c},
       {"epsilon0", PhysConst::ep0},
       {"mu0", PhysConst::mu0},
       {"q_e", PhysConst::q_e},
       {"m_e", PhysConst::m_e},
       {"m_p", PhysConst::m_p},
       {"m_u", PhysConst::m_u},
       {"kb", PhysConst::kb},
       {"pi", MathConst::pi},
      };

    for (auto it = symbols.begin(); it != symbols.end(); ) {
        // Always parsing in double precision avoids potential overflows that may occur when parsing
        // user's expressions because of the limited range of exponentials in single precision
        double v;

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            recursive_symbols.count(*it)==0,
            "Expressions contains recursive symbol "+*it);
        recursive_symbols.insert(*it);
        const bool is_input = queryWithParser(pp_my_constants, it->c_str(), v);
        recursive_symbols.erase(*it);

        if (is_input) {
            parser.setConstant(*it, v);
            it = symbols.erase(it);
            continue;
        }

        auto constant = warpx_constants.find(*it);
        if (constant != warpx_constants.end()) {
          parser.setConstant(*it, constant->second);
          it = symbols.erase(it);
          continue;
        }

        ++it;
    }
    for (auto const& s : symbols) {
        amrex::Abort(Utils::TextMsg::Err("makeParser::Unknown symbol "+s));
    }
    return parser;
}

double
parseStringtoReal(std::string str)
{
    auto parser = makeParser(str, {});
    auto exe = parser.compileHost<0>();
    double result = exe();
    return result;
}

int
parseStringtoInt(std::string str, std::string name)
{
    auto const rval = static_cast<amrex::Real>(parseStringtoReal(str));
    int ival = safeCastToInt(std::round(rval), name);
    return ival;
}
