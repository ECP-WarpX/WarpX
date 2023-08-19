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

    ParmParseWithOptionalGroup(std::string& a_prefix, std::string& a_group);

    bool query (const char* name, bool& ref) const;
    bool query (const char* name, int& ref) const;
    bool query (const char* name, amrex::Real& ref) const;
    bool query (const char* name, std::string& ref) const;

    void get (const char* name, int& ref) const;
    void get (const char* name, long& ref) const;
    void get (const char* name, amrex::Real& ref) const;
    void get (const char* name, std::string& ref) const;
    void get (const char* name, amrex::Vector<amrex::ParticleReal>& ref) const;
    void get (const char* name, amrex::Vector<int>& ref, const int start_ix, const int num_val) const;
    void get (const char* name, amrex::Vector<amrex::Real>& ref, const int start_ix, const int num_val) const;

    void get_long_string (const char* name, std::string& ref) const;
};

#endif
