/* Copyright 2023 David Grote
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARMPARSEWITHOPTIONALGROUP_H_
#define PARMPARSEWITHOPTIONALGROUP_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

struct ParmParseWithOptionalGroup
{
    std::string prefix;
    std::string group;
    std::string prefix_dot_group;

    /* \brief
     * This class wraps around ParmParse and allows querying input parameters
     * with an optional group name.
     * For example, with the base name of "input" and group name of "group",
     * when querying an input name called "variable", this will first check
     * for the input parameter "input.group.variable", and if that is not found,
     * it will check for "input.variable".
     * The group can be an empty string, in which case it will
     * always search for "input.variable".
     *
     * \param[in] a_prefix The base input parameter prefix
     * \param[in] a_group The optional group name in the input parameters
     */
    ParmParseWithOptionalGroup(std::string& a_prefix, std::string& a_group);

    /* Note that this follows the code pattern of ParmParse and writes out
     * each version of the routine rather than using templates */
    bool query (const char* name, bool& ref) const;
    bool query (const char* name, int& ref) const;
    bool query (const char* name, amrex::Real& ref) const;
    bool query (const char* name, std::string& ref) const;

    void get (const char* name, int& ref) const;
    void get (const char* name, long& ref) const;
    void get (const char* name, amrex::Real& ref) const;
    void get (const char* name, std::string& ref) const;
    void get (const char* name, amrex::Vector<amrex::ParticleReal>& ref) const;
    void get (const char* name, amrex::Vector<int>& ref, int start_ix, int num_val) const;
    void get (const char* name, amrex::Vector<amrex::Real>& ref, int start_ix, int num_val) const;

    void get_long_string (const char* name, std::string& ref) const;
};

#endif
