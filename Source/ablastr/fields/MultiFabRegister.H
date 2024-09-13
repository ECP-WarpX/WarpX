/* Copyright 2024 The ABLAST Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Axel Huebl, ...
 */
#ifndef ABLASTR_FIELDS_MF_REGISTER_H
#define ABLASTR_FIELDS_MF_REGISTER_H

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Extension.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <array>
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
     * in separate (i)MultiFab.
     *
     * @todo: synchronize with AMReX "enum class Direction"
     */
    struct Direction
    {
        int dir = 0;

        bool operator<(const Direction& other) const
        {
            return other.dir < this->dir;
        }

        /* TODO: just temporary int compatibility */
        operator int() { return dir; }
    };

    /** A scalar field (a MultiFab)
     *
     * Note: might still have components, e.g., for copies at different times.
     */
    using ScalarField = amrex::MultiFab*;

    /** A read-only scalar field (a MultiFab)
     *
     * Note: might still have components, e.g., for copies at different times.
     */
    using ConstScalarField = amrex::MultiFab const *;

    /** A vector field of three MultiFab
     */
    //using VectorField = ablastr::utils::ConstMap<Direction, amrex::MultiFab *>;
    using VectorField = std::array<amrex::MultiFab *, 3>;

    /** A read-only vector field of three MultiFab
     */
    //using VectorField = ablastr::utils::ConstMap<Direction, amrex::MultiFab const *>;
    using ConstVectorField = std::array<amrex::MultiFab const *, 3>;

    /** A multi-level scalar field
     */
    using MultiLevelScalarField = amrex::Vector<ScalarField>;

    /** A read-only multi-level scalar field
     */
    using ConstMultiLevelScalarField = amrex::Vector<ConstScalarField>;

    /** A multi-level vector field
    */
    using MultiLevelVectorField = amrex::Vector<VectorField>;

    /** A read-only multi-level vector field
    */
    using ConstMultiLevelVectorField = amrex::Vector<ConstVectorField>;

    /** A class to control the lifetime and properties of a MultiFab (field).
     *
     * This class is used to own the lifetime of an amrex::MultiFab and to store
     * associated information around it regarding unique naming, scalar/vector/tensor
     * properties, aliasing, load balancing, etc.
     */
    struct MultiFabOwner
    {
        // TODO: also add iMultiFab via std::variant

        /** owned (i)MultiFab */
        amrex::MultiFab m_mf;

        /** Components (base vector directions) of this (i)MultiFab */
        std::optional<Direction> m_dir = std::nullopt;

        /** the MR level of this (i)MultiFab */
        int m_level = 0;

        /** remake distribution map on load balance, @see amrex::AmrCore::RemakeLevel */
        bool m_remake = true;

        /** redistribute on @see amrex::AmrCore::RemakeLevel */
        bool m_redistribute_on_remake = true;

        /** if m_mf is a non-owning alias, this string tracks the name of the owner */
        std::string m_owner;

        /** Is this part of a vector/tensor? */
        AMREX_INLINE
        bool
        is_vector () const { return m_dir.has_value(); }

        /** Is this an alias MultiFab?
         *
         * If yes, that means we do not own the memory.
         */
        AMREX_INLINE
        bool
        is_alias () const { return !m_owner.empty(); }
    };

    /** This is a register of fields aka amrex::(i)MultiFabs.
     *
     * This is owned by a simulation instance. All used fields should be registered here.
     * Internally, this contains @see MultiFabOwner values.
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
         * @param remake follow the default domain decomposition of the simulation
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
            bool remake = true,
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
            bool remake = true,
            bool redistribute_on_remake = true
        );

        /** Create an alias
         *
         * @param new_name new name
         * @param alias_name owner to alias to
         * @param level ...
         * @param initial_value ...
         * @return
         */
        amrex::MultiFab*
        alias_init (
            std::string new_name,
            std::string alias_name,
            int level,
            std::optional<const amrex::Real> initial_value = std::nullopt
        );

        /** Create an alias
         *
         * @param new_name new name
         * @param alias_name owner to alias to
         * @param dir ...
         * @param level ...
         * @param initial_value ...
         * @return
         */
        amrex::MultiFab*
        alias_init (
                std::string new_name,
                std::string alias_name,
                Direction dir,
                int level,
                std::optional<const amrex::Real> initial_value = std::nullopt
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
        ) const;

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
        ) const;

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
         * body body
         * body
         *
         * @param name ...
         * @param level ...
         * @return ...
         */
        amrex::MultiFab const *
        get (
            std::string name,
            int level
        ) const;

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
        amrex::MultiFab const *
        get (
            std::string name,
            Direction dir,
            int level
        ) const;

        /** title
         *
         * Same as get above, but returns all levels for a name.
         *
         * @param name ...
         * @param finest_level ...
         * @return ...
         */
        MultiLevelScalarField
        get_mr_levels (
            std::string name,
            int finest_level
        );
        ConstMultiLevelScalarField
        get_mr_levels (
            std::string name,
            int finest_level
        ) const;

        /** title
         *
         * Same as get above, but returns all levels for a name.
         *
         * @param name ...
         * @param level ...
         * @return ...
         */
        VectorField
        get_alldirs  (
            std::string name,
            int level
        );
        ConstVectorField
        get_alldirs  (
            std::string name,
            int level
        ) const;

        /** Return a vector field on all MR levels.
         *
         * Out loop: MR levels.
         * Inner loop: directions (components).
         *
         * @param name ...
         * @param finest_level ...
         * @return ...
         */
        MultiLevelVectorField
        get_mr_levels_alldirs (
            std::string name,
            int finest_level
        );
        ConstMultiLevelVectorField
        get_mr_levels_alldirs (
            std::string name,
            int finest_level
        ) const;

        /** title
         *
         * body body
         * body
         *
         * @return ...
         */
        std::vector<std::string>
        list () const;

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

        /** title
         *
         * body body
         * body
         *
         * @param name ...
         * @param dir ...
         * @param level ...
         */
        void
        erase (
            std::string name,
            Direction dir,
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
        ) const;

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
        ) const;

    private:
        amrex::MultiFab *
        internal_get (
                std::string key
        );
        amrex::MultiFab const *
        internal_get (
                std::string key
        ) const;

        /** data storage: ownership and lifetime control */
        std::map<
            std::string,
            MultiFabOwner
        > m_mf_register;
    };

    /** TODO: temporary, remove me
     *
     * @return
     */
    VectorField
    a2m (
        const std::array< std::unique_ptr<amrex::MultiFab>, 3 > & old_vectorfield
    );

    /** TODO: temporary, remove me
     *
     * Convert get_mr_levels_alldirs type to legacy amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >.
     *
     * @return
     */
    MultiLevelVectorField
    va2vm (
        const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& old_vector_on_levels
    );

    /** TODO: temporary, remove me
     *
     * @return
     */
    MultiLevelScalarField
    va2vm (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& old_scalar_on_levels
    );

} // namespace ablastr::fields

#endif  // ABLASTR_FIELDS_MF_REGISTER_H
