/* Copyright 2024 The ABLAST Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Axel Huebl, ...
 */
#include "MultiFabRegister.H"


namespace ablastr::fields
{
    amrex::MultiFab*
    MultiFabRegister::alloc_init (
        std::string name,
        amrex::BoxArray const & ba,
        amrex::DistributionMapping const & dm,
        int ncomp,
        amrex::IntVect const & ngrow,
        int level,
        std::optional<const amrex::Real> initial_value,
        bool redistribute
    )
    {
        name = mf_name(name, level);

        // Checks
        // TODO: does the key already exist? error

        // allocate
        const auto tag = amrex::MFInfo().SetTag(name);
        auto [it, success] = m_mf_register.emplace(
            std::make_pair(
                name,
                MultiFabOwner{{ba, dm, ncomp, ngrow, tag}, redistribute}
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
        std::string /* other_key */,
        int /* other_level */
    )
    {
        // other_key = mf_name(other_key, other_level);

        throw std::runtime_error("MultiFabRegister::alloc_like not yet implemented");

        // Checks
        // TODO: does the key already exist? error
    }

    bool
    MultiFabRegister::contains (
        std::string name,
        int level
    )
    {
        name = mf_name(name, level);

        return m_mf_register.count(name) > 0;
    }

    amrex::MultiFab*
    MultiFabRegister::get (
        std::string name,
        int level
    )
    {
        name = mf_name(name, level);

        if (m_mf_register.count(name) == 0) {
            throw std::runtime_error("MultiFabRegister::get name does not exist in register: " + name);
        }
        amrex::MultiFab & mf = m_mf_register[name].m_mf;

        return &mf;
    }

    std::vector<std::string>
    MultiFabRegister::list ()
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

    std::string
    MultiFabRegister::mf_name (
        std::string name,
        int level
    )
    {
        // Add the suffix "[level=level]"
        return name.append("[level=")
                .append(std::to_string(level))
                .append("]");
    }

} // namespace ablastr::fields
