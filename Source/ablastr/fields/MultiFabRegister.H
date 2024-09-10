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

        bool operator<(const Direction& other) const
        {
            return other.dir < this->dir;
        }
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

        /** Components (base vector directions) of this (i)MultiFab */
        std::optional<Direction> dir = std::nullopt;

        /** the MR level of this (i)MultiFab */
        int level = 0;

        /** redistribute */
        bool redistribute = true;

        /** redistribute on @see amrex::AmrCore::RemakeLevel */
        bool redistribute_on_remake = true;
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
         * @param level ...
         * @param ba ...
         * @param dm ...
         * @param ncomp ...
         * @param ngrow ...
         * @param initial_value ...
         * @param redistribute follow the default domain decomposition of the simulation
         * @param redistribute_on_remake redistribute on @see amrex::AmrCore::RemakeLevel
         * @return pointer to newly allocated MultiFab
         */
        amrex::MultiFab*
        alloc_init (
            std::string name,
            int level,
            amrex::BoxArray const & ba,
            amrex::DistributionMapping const & dm,
            int ncomp,
            amrex::IntVect const & ngrow,
            std::optional<const amrex::Real> initial_value = std::nullopt,
            bool redistribute = true,
            bool redistribute_on_remake = true
        );
        amrex::MultiFab*
        alloc_init (
            std::string name,
            Direction dir,
            int level,
            amrex::BoxArray const & ba,
            amrex::DistributionMapping const & dm,
            int ncomp,
            amrex::IntVect const & ngrow,
            std::optional<const amrex::Real> initial_value = std::nullopt,
            bool redistribute = true,
            bool redistribute_on_remake = true
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
         * @param dir ...
         * @param level ...
         * @return true if contained, otherwise false
         */
        bool
        has (
            std::string name,
            Direction dir,
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

        /** Return a vector field on all MR levels.
         *
         * Out loop: MR levels.
         * Inner loop: directions (components).
         *
         * @param name ...
         * @param finest_level ...
         * @return ...
         */
        std::vector<
            std::map<
                Direction,
                amrex::MultiFab*
            >
        >
        get_mr_levels_alldirs (
            std::string name,
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

        /** Remake all (i)MultiFab with a new distribution mapping.
         *
         * If redistribute is true, we also copy from the old data into the new.
         *
         * @param level the MR level to erase all MultiFabs from
         * @param new_dm new distribution mapping
         */
        void
        remake_level (
            int other_level,
            amrex::DistributionMapping const & new_dm
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

    /** TODO: temporary, remove me
     *
     * Convert get_mr_levels_alldirs type to legacy amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >.
     *
     * @return
     */
    std::vector<
        std::map<
                Direction,
                amrex::MultiFab*
        >
    >
    va2vm (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > const & old_vector_on_levels
    );

} // namespace ablastr::fields

#endif  // ABLASTR_FIELDS_MF_REGISTER_H
