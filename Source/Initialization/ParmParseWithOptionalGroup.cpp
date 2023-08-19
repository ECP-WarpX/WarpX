/* Copyright 2023 David Grote
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ParmParseWithOptionalGroup.H"
#include "Utils/Parser/ParserUtils.H"

#include <AMReX_ParmParse.H>

#include <vector>

ParmParseWithOptionalGroup::ParmParseWithOptionalGroup(std::string& a_prefix, std::string& a_group):
    prefix{a_prefix}, group{a_group}
{
    if (group.empty()) {
        prefix_dot_group = prefix;
    } else {
        prefix_dot_group = prefix + "." + group;
    }
}

bool ParmParseWithOptionalGroup::query (const char* name, bool& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = utils::parser::queryWithParser(pp_with_group, name, ref);
    if (!is_specified && !group.empty()) {
        // If it wasn't found, query it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        is_specified = utils::parser::queryWithParser(pp_prefix, name, ref);
    }
    return is_specified;
}

bool ParmParseWithOptionalGroup::query (const char* name, int& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = utils::parser::queryWithParser(pp_with_group, name, ref);
    if (!is_specified && !group.empty()) {
        // If it wasn't found, query it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        is_specified = utils::parser::queryWithParser(pp_prefix, name, ref);
    }
    return is_specified;
}

bool ParmParseWithOptionalGroup::query (const char* name, amrex::Real& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = utils::parser::queryWithParser(pp_with_group, name, ref);
    if (!is_specified && !group.empty()) {
        // If it wasn't found, query it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        is_specified = utils::parser::queryWithParser(pp_prefix, name, ref);
    }
    return is_specified;
}

bool ParmParseWithOptionalGroup::query (const char* name, std::string& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = pp_with_group.query(name, ref);
    if (!is_specified && !group.empty()) {
        // If it wasn't found, query it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        is_specified = pp_prefix.query(name, ref);
    }
    return is_specified;
}

void ParmParseWithOptionalGroup::get_long_string (const char* name, std::string& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = pp_with_group.contains(name);
    if (is_specified) {
        utils::parser::Store_parserString(pp_with_group, name, ref);
    } else {
        // If it wasn't found, query it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        utils::parser::Store_parserString(pp_prefix, name, ref);
    }

}

void ParmParseWithOptionalGroup::get (const char* name, int& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = utils::parser::queryWithParser(pp_with_group, name, ref);
    if (!is_specified) {
        // If it wasn't found, get it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        utils::parser::getWithParser(pp_prefix, name, ref);
    }
}

void ParmParseWithOptionalGroup::get (const char* name, long& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = utils::parser::queryWithParser(pp_with_group, name, ref);
    if (!is_specified) {
        // If it wasn't found, get it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        utils::parser::getWithParser(pp_prefix, name, ref);
    }
}

void ParmParseWithOptionalGroup::get (const char* name, amrex::Real& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = utils::parser::queryWithParser(pp_with_group, name, ref);
    if (!is_specified) {
        // If it wasn't found, get it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        utils::parser::getWithParser(pp_prefix, name, ref);
    }
}

void ParmParseWithOptionalGroup::get (const char* name, std::string& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = pp_with_group.query(name, ref);
    if (!is_specified) {
        // If it wasn't found, get it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        pp_prefix.get(name, ref);
    }
}

void ParmParseWithOptionalGroup::get (const char* name, amrex::Vector<amrex::ParticleReal>& ref) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_specified = utils::parser::queryArrWithParser(pp_with_group, name, ref);
    if (!is_specified) {
        // If it wasn't found, get it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        utils::parser::getArrWithParser(pp_prefix, name, ref);
    }
}

void ParmParseWithOptionalGroup::get (const char* name, amrex::Vector<int>& ref, const int start_ix, const int num_val) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_with_suffix = utils::parser::queryArrWithParser(pp_with_group, name, ref, start_ix, num_val);
    if (!is_with_suffix) {
        // If it wasn't found, get it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        utils::parser::getArrWithParser(pp_prefix, name, ref, start_ix, num_val);
    }
}

void ParmParseWithOptionalGroup::get (const char* name, amrex::Vector<amrex::Real>& ref, const int start_ix, const int num_val) const
{
    // First, query name with the group name
    const amrex::ParmParse pp_with_group(prefix_dot_group);
    bool is_with_suffix = utils::parser::queryArrWithParser(pp_with_group, name, ref, start_ix, num_val);
    if (!is_with_suffix) {
        // If it wasn't found, get it with only the prefix
        const amrex::ParmParse pp_prefix(prefix);
        utils::parser::getArrWithParser(pp_prefix, name, ref, start_ix, num_val);
    }
}
