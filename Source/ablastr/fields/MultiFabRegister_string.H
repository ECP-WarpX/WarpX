/* Copyright 2024 The ABLASTR Community
 *
 * This file is part of ABLASTR.
 * It adds unique string key access to the MultiFabRegister.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Axel Huebl
 */
#ifndef ABLASTR_FIELDS_MF_REGISTER_STRING_H
#define ABLASTR_FIELDS_MF_REGISTER_STRING_H

#include "MultiFabRegister.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <optional>
#include <string>


namespace ablastr::fields
{
    template<>
    amrex::MultiFab*
    MultiFabRegister::alloc_init<std::string const &> (
        std::string const & name,
        int level,
        amrex::BoxArray const & ba,
        amrex::DistributionMapping const & dm,
        int ncomp,
        amrex::IntVect const & ngrow,
        std::optional<amrex::Real const> initial_value,
        bool remake,
        bool redistribute_on_remake
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::alloc_init<std::string const &> (
        std::string const & name,
        Direction dir,
        int level,
        amrex::BoxArray const & ba,
        amrex::DistributionMapping const & dm,
        int ncomp,
        amrex::IntVect const & ngrow,
        std::optional<amrex::Real const> initial_value,
        bool remake,
        bool redistribute_on_remake
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::alias_init<std::string const &, std::string const &> (
        std::string const & new_name,
        std::string const & alias_name,
        int level,
        std::optional<amrex::Real const> initial_value
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::alias_init<std::string const &, std::string const &> (
        std::string const & new_name,
        std::string const & alias_name,
        Direction dir,
        int level,
        std::optional<amrex::Real const> initial_value
    );

    template<>
    bool
    MultiFabRegister::has<std::string const &> (
        std::string const & name,
        int level
    ) const;

    template<>
    bool
    MultiFabRegister::has<std::string const &> (
        std::string const & name,
        Direction dir,
        int level
    ) const;

    template<>
    bool
    MultiFabRegister::has_vector<std::string const &> (
        std::string const & name,
        int level
    ) const;

    template<>
    amrex::MultiFab*
    MultiFabRegister::get<std::string const &> (
        std::string const & name,
        int level
    );

    template<>
    amrex::MultiFab*
    MultiFabRegister::get<std::string const &> (
        std::string const & name,
        Direction dir,
        int level
    );

    template<>
    amrex::MultiFab const*
    MultiFabRegister::get<std::string const &> (
        std::string const & name,
        int level
    ) const;

    template<>
    amrex::MultiFab const *
    MultiFabRegister::get<std::string const &> (
        std::string const & name,
        Direction dir,
        int level
    ) const;

    template<>
    MultiLevelScalarField
    MultiFabRegister::get_mr_levels<std::string const &> (
        std::string const & name,
        int finest_level
    );

    template<>
    ConstMultiLevelScalarField
    MultiFabRegister::get_mr_levels<std::string const &> (
        std::string const & name,
        int finest_level
    ) const;

    template<>
    VectorField
    MultiFabRegister::get_alldirs<std::string const &> (
        std::string const & name,
        int level
    );

    template<>
    ConstVectorField
    MultiFabRegister::get_alldirs<std::string const &> (
        std::string const & name,
        int level
    ) const;

    template<>
    MultiLevelVectorField
    MultiFabRegister::get_mr_levels_alldirs<std::string const &> (
        std::string const & name,
        int finest_level
    );

    template<>
    ConstMultiLevelVectorField
    MultiFabRegister::get_mr_levels_alldirs<std::string const &> (
        std::string const & name,
        int finest_level
    ) const;

    template<>
    void
    MultiFabRegister::erase<std::string const &> (
        std::string const & name,
        int level
    );

    template<>
    void
    MultiFabRegister::erase<std::string const &> (
        std::string const & name,
        Direction dir,
        int level
    );

} // namespace ablastr::fields

#endif  // ABLASTR_FIELDS_MF_REGISTER_STRING_H
