/* Copyright 2024 The ABLAST Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Axel Huebl, ...
 */
#ifndef ABLASTR_FIELDS_MF_REGISTER_H
#define ABLASTR_FIELDS_MF_REGISTER_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <stdexcept>
#include <vector>


namespace ablastr::fields
{
    struct MultiFabOwner
    {
        // TODO: also add iMultiFab via std::variant

        /** owned (i)MultiFab data */
        amrex::MultiFab m_mf;

        /** redistribute */
        bool redistribute = true;
    };

    /** This is a register of fields aka amrex::(i)MultiFabs.
     *
     * This is owned by a simulation instance. All used fields should be registered here.
     */
    struct MultiFabRegister
    {
        // Avoid accidental copies when passing to member functions
        MultiFabRegister() = default;
        MultiFabRegister(MultiFabRegister const &) = delete;
        MultiFabRegister(MultiFabRegister&&) = delete;
        MultiFabRegister& operator=(MultiFabRegister const &) = delete;
        MultiFabRegister& operator=(MultiFabRegister&&) = delete;
        ~MultiFabRegister() = default;

        /** title
         *
         * body body
         * body
         * body
         *
         * @param name ...
         * @param ba ...
         * @param dm ...
         * @param ncomp ...
         * @param ngrow ...
         * @param level ...
         * @param initial_value ...
         * @param redistribute follow the default domain decomposition of the simulation
         * @return pointer to newly allocated MultiFab
         */
        amrex::MultiFab*
        alloc_init (
            std::string name,
            amrex::BoxArray const & ba,
            amrex::DistributionMapping const & dm,
            int ncomp,
            amrex::IntVect const & ngrow,
            int level,
            std::optional<const amrex::Real> initial_value = std::nullopt,
            bool redistribute = true
        );

        /** title
         *
         * body body
         * body
         *
         * @param name ...
         * @param other_level ...
         */
        void
        alloc_like (
            std::string other_key,
            int other_level
        );

        /** title
         *
         * body body
         * body
         *
         * @param name ...
         * @param level ...
         * @return true if contained, otherwise false
         */
        bool
        contains (
            std::string name,
            int level
        );

        /** title
         *
         * body body
         * body
         *
         * @param name ...
         * @param level ...
         * @return ...
         */
        amrex::MultiFab*
        get (
            std::string name,
            int level
        );

        /** title
         *
         * body body
         * body
         *
         * @return ...
         */
        std::vector<std::string>
        list ();

        /** title
         *
         * body body
         * body
         *
         * @param name ...
         * @param level ...
         * @return ...
         */
        void
        erase (
            std::string name,
            int level
        );

        /** title
         *
         * body body
         * body
         *
         * @param name ...
         * @param level ...
         * @return ...
         */
        std::string
        mf_name (
            std::string name,
            int level
        );

    private:
        /** data storage: ownership and lifetime control */
        std::map<
            std::string,
            MultiFabOwner
        > m_mf_register;
    };

} // namespace ablastr::fields

#endif  // ABLASTR_FIELDS_MF_REGISTER_H
