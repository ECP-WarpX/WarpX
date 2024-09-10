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

    /** Components (base vector directions) of vector/tensor fields.
     *
     * Because of different staggering, the components of vector/tensor fields are stored
     * in separate (i)MultiFab. This
     */
    struct Direction
    {
        int dir = 0;
    };

    /** title
     *
     * body body
     * body
     * body
     */
    struct MultiFabOwner
    {
        // TODO: also add iMultiFab via std::variant

        /** owned (i)MultiFab data */
        amrex::MultiFab m_mf;

        /** the MR level of this (i)MultiFab */
        int level = 0;

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
         * @param other_name ...
         * @param other_level ...
         */
        void
        alloc_like (
            std::string other_name,
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
        has (
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
         * @param name ...
         * @param dir ...
         * @param level ...
         * @return ...
         */
        amrex::MultiFab*
        get (
            std::string name,
            Direction dir,
            int level
        );

        /** title
         *
         * Same as get above, but returns all levels for a name.
         *
         * @param name ...
         * @param finest_level ...
         * @return ...
         */
        std::vector<amrex::MultiFab*>
        get_mr_levels (
            std::string name,
            int finest_level
        );

        /** title
         *
         * Same as get above, but returns all levels for a name.
         *
         * @param name ...
         * @param dir ...
         * @param finest_level ...
         * @return ...
         */
        std::vector<amrex::MultiFab*>
        get_mr_levels (
            std::string name,
            Direction dir,
            int finest_level
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
         */
        void
        erase (
            std::string name,
            int level
        );

        /** Erase all MultiFabs on a specific MR level.
         *
         * Calls @see erase for all MultiFabs on a specific level.
         *
         * @param level the MR level to erase all MultiFabs from
         */
        void
        clear_level (
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

        /** title
         *
         * body body
         * body
         *
         * @param name ...
         * @param dir ...
         * @param level ...
         * @return ...
         */
        std::string
        mf_name (
            std::string name,
            Direction dir,
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
