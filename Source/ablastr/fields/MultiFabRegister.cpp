/* Copyright 2024 The ABLAST Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Axel Huebl, ...
 */
#include "MultiFabRegister.H"

#include "ablastr/utils/TextMsg.H"

#include <algorithm>


namespace ablastr::fields
{
    amrex::MultiFab*
    MultiFabRegister::alloc_init (
        std::string name,
        int level,
        amrex::BoxArray const & ba,
        amrex::DistributionMapping const & dm,
        int ncomp,
        amrex::IntVect const & ngrow,
        std::optional<const amrex::Real> initial_value,
        bool remake,
        bool redistribute_on_remake
    )
    {
        // checks
        if (has(name, level)) {
            throw std::runtime_error("MultiFabRegister::alloc_init failed because " + name + " already exists.");
        }

        // fully qualified name
        name = mf_name(name, level);

        // allocate
        const auto tag = amrex::MFInfo().SetTag(name);
        auto [it, success] = m_mf_register.emplace(
            std::make_pair(
                name,
                MultiFabOwner{
                    {ba, dm, ncomp, ngrow, tag},
                    std::nullopt,  // scalar: no direction
                    level,
                    remake,
                    redistribute_on_remake,
                    ""   // we own the memory
                }
            )
        );
        if (!success) {
            throw std::runtime_error("MultiFabRegister::alloc_init failed for " + name);
        }

        // a short-hand alias for the code below
        amrex::MultiFab & mf = it->second.m_mf;

        // initialize with value
        if (initial_value) {
            mf.setVal(*initial_value);
        }

        return &mf;
    }

    amrex::MultiFab*
    MultiFabRegister::alloc_init (
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
    )
    {
        // checks
        if (has(name, dir, level)) {
            throw std::runtime_error(
                "MultiFabRegister::alloc_init failed because " +
                mf_name(name, dir, level) +
                " already exists."
            );
        }

        // fully qualified name
        name = mf_name(name, dir, level);

        // allocate
        const auto tag = amrex::MFInfo().SetTag(name);
        auto [it, success] = m_mf_register.emplace(
            std::make_pair(
                name,
                MultiFabOwner{
                    {ba, dm, ncomp, ngrow, tag},
                    dir,
                    level,
                    remake,
                    redistribute_on_remake,
                    ""   // we own the memory
                }
            )
        );
        if (!success) {
            throw std::runtime_error("MultiFabRegister::alloc_init failed for " + name);
        }

        // a short-hand alias for the code below
        amrex::MultiFab & mf = it->second.m_mf;

        // initialize with value
        if (initial_value) {
            mf.setVal(*initial_value);
        }

        return &mf;
    }

    void
    MultiFabRegister::alloc_like (
        std::string /* other_name */,
        int /* other_level */
    )
    {
        // other_name = mf_name(other_name, other_level);

        throw std::runtime_error("MultiFabRegister::alloc_like not yet implemented");

        // Checks
        // TODO: does the other_name already exist? error
    }

    amrex::MultiFab*
    MultiFabRegister::alias_init (
        std::string new_name,
        std::string alias_name,
        int level,
        std::optional<const amrex::Real> initial_value
    )
    {
        // checks
        if (has(new_name, level)) {
            throw std::runtime_error(
                "MultiFabRegister::alias_init failed because " +
                mf_name(new_name, level) +
                " already exists."
            );
        }
        if (!has(alias_name, level)) {
            throw std::runtime_error(
                    "MultiFabRegister::alias_init failed because " +
                    mf_name(alias_name, level) +
                    " does not exist."
            );
        }

        // fully qualified name
        new_name = mf_name(new_name, level);
        alias_name = mf_name(alias_name, level);

        MultiFabOwner & alias = m_mf_register[alias_name];
        amrex::MultiFab & mf_alias = alias.m_mf;

        // allocate
        auto [it, success] = m_mf_register.emplace(
            std::make_pair(
                new_name,
                MultiFabOwner{
                    {mf_alias, amrex::make_alias, 0, mf_alias.nComp()},
                    std::nullopt,  // scalar: no direction
                    level,
                    alias.m_remake,
                    alias.m_redistribute_on_remake,
                    alias_name
                }
            )
        );
        if (!success) {
            throw std::runtime_error("MultiFabRegister::alias_init failed for " + new_name);
        }

        // a short-hand alias for the code below
        amrex::MultiFab & mf = it->second.m_mf;

        // initialize with value
        if (initial_value) {
            mf.setVal(*initial_value);
        }

        return &mf;
    }

    amrex::MultiFab*
    MultiFabRegister::alias_init (
            std::string new_name,
            std::string alias_name,
            Direction dir,
            int level,
            std::optional<const amrex::Real> initial_value
    )
    {
        // checks
        if (has(new_name, dir, level)) {
            throw std::runtime_error(
                "MultiFabRegister::alias_init failed because " +
                mf_name(new_name, dir, level) +
                " already exists."
            );
        }
        if (!has(alias_name, dir, level)) {
            throw std::runtime_error(
                "MultiFabRegister::alias_init failed because " +
                mf_name(alias_name, dir, level) +
                " does not exist."
            );
        }

        // fully qualified name
        new_name = mf_name(new_name, dir, level);
        alias_name = mf_name(alias_name, dir, level);

        MultiFabOwner & alias = m_mf_register[alias_name];
        amrex::MultiFab & mf_alias = alias.m_mf;

        // allocate
        auto [it, success] = m_mf_register.emplace(
            std::make_pair(
                new_name,
                MultiFabOwner{
                    {mf_alias, amrex::make_alias, 0, mf_alias.nComp()},
                    dir,
                    level,
                    alias.m_remake,
                    alias.m_redistribute_on_remake,
                    alias_name
                }
            )
        );
        if (!success) {
            throw std::runtime_error("MultiFabRegister::alias_init failed for " + new_name);
        }

        // a short-hand alias for the code below
        amrex::MultiFab & mf = it->second.m_mf;

        // initialize with value
        if (initial_value) {
            mf.setVal(*initial_value);
        }

        return &mf;
    }

    void
    MultiFabRegister::remake_level (
        int level,
        amrex::DistributionMapping const & new_dm
    )
    {
        // Owning MultiFabs
        for (auto & element : m_mf_register )
        {
            MultiFabOwner & mf_owner = element.second;

            // keep distribution map as it is?
            if (!mf_owner.m_remake) {
                continue;
            }

            // remake MultiFab with new distribution map
            if (mf_owner.m_level == level && !mf_owner.is_alias()) {
                amrex::MultiFab & mf = mf_owner.m_mf;
                amrex::IntVect const & ng = mf.nGrowVect();
                const auto tag = amrex::MFInfo().SetTag(mf.tags()[0]);
                amrex::MultiFab new_mf(mf.boxArray(), new_dm, mf.nComp(), ng, tag);

                // copy data to new MultiFab: Only done for persistent data like E and B field, not for
                // temporary things like currents, etc.
                if (mf_owner.m_redistribute_on_remake) {
                    new_mf.Redistribute(mf, 0, 0, mf.nComp(), ng);
                }

                // replace old MultiFab with new one, deallocate old one
                mf_owner.m_mf = std::move(new_mf);
            }
        }

        // Aliases
        for (auto & element : m_mf_register )
        {
            MultiFabOwner & mf_owner = element.second;

            // keep distribution map as it is?
            if (!mf_owner.m_remake) {
                continue;
            }

            if (mf_owner.m_level == level && mf_owner.is_alias()) {
                amrex::MultiFab & mf = m_mf_register[mf_owner.m_owner].m_mf;
                amrex::MultiFab new_mf(mf, amrex::make_alias, 0, mf.nComp());

                // no copy via Redistribute: the owner was already redistributed

                // replace old MultiFab with new one, deallocate old one
                mf_owner.m_mf = std::move(new_mf);
            }
        }
    }

    bool
    MultiFabRegister::has (
        std::string name,
        int level
    ) const
    {
        name = mf_name(name, level);

        return m_mf_register.count(name) > 0;
    }

    bool
    MultiFabRegister::has (
        std::string name,
        Direction dir,
        int level
    ) const
    {
        name = mf_name(name, dir, level);

        return m_mf_register.count(name) > 0;
    }

    amrex::MultiFab*
    MultiFabRegister::internal_get (
        std::string key
    )
    {
        if (m_mf_register.count(key) == 0) {
            // FIXME: temporary, throw a std::runtime_error
            // throw std::runtime_error("MultiFabRegister::get name does not exist in register: " + key);
            return nullptr;
        }
        amrex::MultiFab & mf = m_mf_register.at(key).m_mf;

        return &mf;
    }

    amrex::MultiFab const *
    MultiFabRegister::internal_get (
        std::string key
    ) const
    {
        if (m_mf_register.count(key) == 0) {
            // FIXME: temporary, throw a std::runtime_error
            // throw std::runtime_error("MultiFabRegister::get name does not exist in register: " + key);
            return nullptr;
        }
        amrex::MultiFab const & mf = m_mf_register.at(key).m_mf;

        return &mf;
    }

    amrex::MultiFab*
    MultiFabRegister::get (
        std::string name,
        int level
    )
    {
        name = mf_name(name, level);
        return internal_get(name);
    }

    amrex::MultiFab*
    MultiFabRegister::get (
        std::string name,
        Direction dir,
        int level
    )
    {
        name = mf_name(name, dir, level);
        return internal_get(name);
    }

    amrex::MultiFab const *
    MultiFabRegister::get (
        std::string name,
        int level
    ) const
    {
        name = mf_name(name, level);
        return internal_get(name);
    }

    amrex::MultiFab const *
    MultiFabRegister::get (
        std::string name,
        Direction dir,
        int level
    ) const
    {
        name = mf_name(name, dir, level);
        return internal_get(name);
    }

    MultiLevelScalarField
    MultiFabRegister::get_mr_levels (
        std::string name,
        int finest_level
    )
    {
        MultiLevelScalarField field_on_level;
        field_on_level.reserve(finest_level+1);
        for (int lvl = 0; lvl <= finest_level; lvl++)
        {
            field_on_level.push_back(get(name, lvl));
        }
        return field_on_level;
    }

    ConstMultiLevelScalarField
    MultiFabRegister::get_mr_levels (
        std::string name,
        int finest_level
    ) const
    {
        ConstMultiLevelScalarField field_on_level;
        field_on_level.reserve(finest_level+1);
        for (int lvl = 0; lvl <= finest_level; lvl++)
        {
            field_on_level.push_back(get(name, lvl));
        }
        return field_on_level;
    }

    VectorField
    MultiFabRegister::get_alldirs  (
        std::string name,
        int level
    )
    {
        // TODO: Technically, we should search field_on_level via std::unique_copy
        std::vector<Direction> all_dirs = {Direction{0}, Direction{1}, Direction{2}};

        // insert a new level
        VectorField vectorField;

        // insert components
        for (Direction dir : all_dirs)
        {
            vectorField[dir] = get(name, dir, level);
        }
        return vectorField;
    }

    ConstVectorField
    MultiFabRegister::get_alldirs  (
        std::string name,
        int level
    ) const
    {
        // TODO: Technically, we should search field_on_level via std::unique_copy
        std::vector<Direction> all_dirs = {Direction{0}, Direction{1}, Direction{2}};

        // insert a new level
        ConstVectorField vectorField;

        // insert components
        for (Direction dir : all_dirs)
        {
            vectorField[dir] = get(name, dir, level);
        }
        return vectorField;
    }

    MultiLevelVectorField
    MultiFabRegister::get_mr_levels_alldirs  (
        std::string name,
        int finest_level
    )
    {
        MultiLevelVectorField field_on_level;
        field_on_level.reserve(finest_level+1);

        // TODO: Technically, we should search field_on_level via std::unique_copy
        std::vector<Direction> all_dirs = {Direction{0}, Direction{1}, Direction{2}};

        for (int lvl = 0; lvl <= finest_level; lvl++)
        {
            // insert a new level
            field_on_level.push_back(VectorField{});

            // insert components
            for (Direction dir : all_dirs)
            {
                field_on_level[lvl][dir] = get(name, dir, lvl);
            }
        }
        return field_on_level;
    }

    ConstMultiLevelVectorField
    MultiFabRegister::get_mr_levels_alldirs  (
        std::string name,
        int finest_level
    ) const
    {
        ConstMultiLevelVectorField field_on_level;
        field_on_level.reserve(finest_level+1);

        // TODO: Technically, we should search field_on_level via std::unique_copy
        std::vector<Direction> all_dirs = {Direction{0}, Direction{1}, Direction{2}};

        for (int lvl = 0; lvl <= finest_level; lvl++)
        {
            // insert a new level
            field_on_level.push_back(ConstVectorField{});

            // insert components
            for (Direction dir : all_dirs)
            {
                field_on_level[lvl][dir] = get(name, dir, lvl);
            }
        }
        return field_on_level;
    }

    std::vector<std::string>
    MultiFabRegister::list () const
    {
        std::vector<std::string> names;
        names.reserve(m_mf_register.size());
        for (auto const & str : m_mf_register) { names.push_back(str.first); }

        return names;
    }

    void
    MultiFabRegister::erase (
        std::string name,
        int level
    )
    {
        name = mf_name(name, level);

        if (m_mf_register.count(name) != 1) {
            throw std::runtime_error("MultiFabRegister::remove name does not exist in register: " + name);
        }
        m_mf_register.erase(name);
    }

    void
    MultiFabRegister::erase (
        std::string name,
        Direction dir,
        int level
    )
    {
        name = mf_name(name, dir, level);

        if (m_mf_register.count(name) != 1) {
            throw std::runtime_error("MultiFabRegister::remove name does not exist in register: " + name);
        }
        m_mf_register.erase(name);
    }

    void
    MultiFabRegister::clear_level (
        int level
    )
    {
        // C++20: Replace with std::erase_if
        for (auto first = m_mf_register.begin(), last = m_mf_register.end(); first != last;)
        {
            if (first->second.m_level == level)
                first = m_mf_register.erase(first);
            else
                ++first;
        }
    }

    std::string
    MultiFabRegister::mf_name (
        std::string name,
        int level
    ) const
    {
        // Add the suffix "[level=level]"
        return name.append("[level=")
                .append(std::to_string(level))
                .append("]");
    }

    std::string
    MultiFabRegister::mf_name (
        std::string name,
        Direction dir,
        int level
    ) const
    {
        // Add the suffix "[dir=dir]"
        return mf_name(
            name
            .append("[dir=")
            .append(std::to_string(dir.dir))
            .append("]"),
            level
        );
    }

    VectorField
    a2m (
        const std::array< std::unique_ptr<amrex::MultiFab>, 3 > & old_vectorfield
    )
    {
        std::vector<Direction> all_dirs = {Direction{0}, Direction{1}, Direction{2}};

        VectorField field_on_level;

        // insert components
        for (auto dir : {0, 1, 2})
        {
            field_on_level[Direction{dir}] = old_vectorfield[dir].get();
        }
        return field_on_level;
    }

    MultiLevelVectorField
    va2vm (
        const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& old_vector_on_levels
    )
    {
        int const finest_level = old_vector_on_levels.size() - 1u;

        MultiLevelVectorField field_on_level;
        field_on_level.reserve(finest_level+1);

        std::vector<Direction> all_dirs = {Direction{0}, Direction{1}, Direction{2}};

        for (int lvl = 0; lvl <= finest_level; lvl++)
        {
            // insert a new level
            field_on_level.push_back(VectorField{});

            // insert components
            for (auto dir : {0, 1, 2})
            {
                field_on_level[lvl][Direction{dir}] = old_vector_on_levels[lvl][dir].get();
            }
        }
        return field_on_level;
    }

    MultiLevelScalarField
    va2vm (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& old_scalar_on_levels
    )
    {
        int const finest_level = old_scalar_on_levels.size() - 1u;

        MultiLevelScalarField field_on_level;
        field_on_level.reserve(finest_level+1);

        for (int lvl = 0; lvl <= finest_level; lvl++)
        {
            // insert a scalar field on a level
            field_on_level.push_back(ScalarField{old_scalar_on_levels[lvl].get()});

        }
        return field_on_level;
    }
} // namespace ablastr::fields
