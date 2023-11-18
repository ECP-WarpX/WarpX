/* Copyright 2023 David Grote
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ParmParseWithOptionalGroup.H"
#include "Utils/Parser/ParserUtils.H"

#include <vector>

ParmParseWithOptionalGroup::ParmParseWithOptionalGroup(std::string const& a_prefix, std::string const& a_group):
    prefix{a_prefix}, group{a_group}
{
    if (group.empty()) {
        prefix_dot_group = prefix;
    } else {
        prefix_dot_group = prefix + "." + group;
    }
}

const amrex::ParmParse ParmParseWithOptionalGroup::SelectParmParse(const char *name) const
{
    // Check if the name is found with and without the group.
    const amrex::ParmParse pp_prefix(prefix);
    const bool is_specified_without_group = pp_prefix.contains(name);

    const amrex::ParmParse pp_with_group(prefix_dot_group);
    const bool is_specified_with_group = (group.empty() ? false : pp_with_group.contains(name));

    if (is_specified_without_group && !is_specified_with_group) {
        // If found without the group but not with the group, then use the one without the group.
        return pp_prefix;
    } else {
        // Otherwise, use the one with the group even if not found, in which case an exception may be raised.
        return pp_with_group;
    }
}

bool ParmParseWithOptionalGroup::queryWithParser (const char* name, bool& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    return utils::parser::queryWithParser(pp, name, ref);
}

bool ParmParseWithOptionalGroup::queryWithParser (const char* name, int& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    return utils::parser::queryWithParser(pp, name, ref);
}

bool ParmParseWithOptionalGroup::queryWithParser (const char* name, double& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    return utils::parser::queryWithParser(pp, name, ref);
}

bool ParmParseWithOptionalGroup::queryWithParser (const char* name, float& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    return utils::parser::queryWithParser(pp, name, ref);
}

bool ParmParseWithOptionalGroup::query (const char* name, std::string& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    return pp.query(name, ref);
}

void ParmParseWithOptionalGroup::get_long_string (const char* name, std::string& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::Store_parserString(pp, name, ref);
}

void ParmParseWithOptionalGroup::getWithParser (const char* name, int& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getWithParser(pp, name, ref);
}

void ParmParseWithOptionalGroup::getWithParser (const char* name, long& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getWithParser(pp, name, ref);
}

void ParmParseWithOptionalGroup::getWithParser (const char* name, double& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getWithParser(pp, name, ref);
}

void ParmParseWithOptionalGroup::getWithParser (const char* name, float& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getWithParser(pp, name, ref);
}

void ParmParseWithOptionalGroup::get (const char* name, std::string& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    pp.get(name, ref);
}

void ParmParseWithOptionalGroup::getArrWithParser (const char* name, amrex::Vector<double>& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getArrWithParser(pp, name, ref);
}

void ParmParseWithOptionalGroup::getArrWithParser (const char* name, amrex::Vector<float>& ref) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getArrWithParser(pp, name, ref);
}

void ParmParseWithOptionalGroup::getArrWithParser (const char* name, amrex::Vector<int>& ref, int start_ix, int num_val) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getArrWithParser(pp, name, ref, start_ix, num_val);
}

void ParmParseWithOptionalGroup::getArrWithParser (const char* name, amrex::Vector<double>& ref, int start_ix, int num_val) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getArrWithParser(pp, name, ref, start_ix, num_val);
}

void ParmParseWithOptionalGroup::getArrWithParser (const char* name, amrex::Vector<float>& ref, int start_ix, int num_val) const
{
    const amrex::ParmParse pp = SelectParmParse(name);
    utils::parser::getArrWithParser(pp, name, ref, start_ix, num_val);
}
