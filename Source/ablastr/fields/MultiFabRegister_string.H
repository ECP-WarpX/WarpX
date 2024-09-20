/* Copyright 2024 The ABLAST Community
 *
 * This file is part of WarpX.
 * It adds unique string key access to the MultiFabRegister.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Axel Huebl
 */
#ifndef ABLASTR_FIELDS_MF_REGISTER_STRING_H
#define ABLASTR_FIELDS_MF_REGISTER_STRING_H

#include "MultiFabRegister.H"


namespace ablastr::fields
{
    template<>
    amrex::MultiFab*
    MultiFabRegister::alloc_init<std::string> (
        std::string name,
        int level,
        amrex::BoxArray const & ba,
        amrex::DistributionMapping const & dm,
        int ncomp,
        amrex::IntVect const & ngrow,
        std::optional<const amrex::Real> initial_value,
        bool remake,
        bool redistribute_on_remake
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::alloc_init<std::string> (
        std::string name,
        Direction dir,
        int level,
        amrex::BoxArray const & ba,
        amrex::DistributionMapping const & dm,
        int ncomp,
        amrex::IntVect const & ngrow,
        std::optional<const amrex::Real> initial_value,
        bool remake,
        bool redistribute_on_remake
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::alias_init<std::string, std::string> (
        std::string new_name,
        std::string alias_name,
        int level,
        std::optional<const amrex::Real> initial_value
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::alias_init<std::string, std::string> (
        std::string new_name,
        std::string alias_name,
        Direction dir,
        int level,
        std::optional<const amrex::Real> initial_value
    );

    template<>
    bool
    MultiFabRegister::has<std::string> (
        std::string name,
        int level
    ) const;

    template<>
    bool
    MultiFabRegister::has<std::string> (
        std::string name,
        Direction dir,
        int level
    ) const;

    template<>
    bool
    MultiFabRegister::has_vector<std::string> (
        std::string name,
        int level
    ) const;

    template<>
    amrex::MultiFab*
    MultiFabRegister::get<std::string> (
        std::string name,
        int level
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::get<std::string> (
        std::string name,
        Direction dir,
        int level
    );

    template<>
    amrex::MultiFab const*
    MultiFabRegister::get<std::string> (
        std::string name,
        int level
    ) const;

    template<>
    amrex::MultiFab const *
    MultiFabRegister::get<std::string> (
        std::string name,
        Direction dir,
        int level
    ) const;

    template<>
    MultiLevelScalarField
    MultiFabRegister::get_mr_levels<std::string> (
        std::string name,
        int finest_level
    );

    template<>
    ConstMultiLevelScalarField
    MultiFabRegister::get_mr_levels<std::string> (
        std::string name,
        int finest_level
    ) const;

    template<>
    VectorField
    MultiFabRegister::get_alldirs<std::string> (
        std::string name,
        int level
    );

    template<>
    ConstVectorField
    MultiFabRegister::get_alldirs<std::string> (
        std::string name,
        int level
    ) const;

    template<>
    MultiLevelVectorField
    MultiFabRegister::get_mr_levels_alldirs<std::string> (
        std::string name,
        int finest_level
    );

    template<>
    ConstMultiLevelVectorField
    MultiFabRegister::get_mr_levels_alldirs<std::string> (
        std::string name,
        int finest_level
    ) const;

    template<>
    void
    MultiFabRegister::erase<std::string> (
        std::string name,
        int level
    );

    template<>
    void
    MultiFabRegister::erase<std::string> (
        std::string name,
        Direction dir,
        int level
    );

} // namespace ablastr::fields

#endif  // ABLASTR_FIELDS_MF_REGISTER_STRING_H
