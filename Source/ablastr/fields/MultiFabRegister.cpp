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
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        const int ncomp,
        const amrex::IntVect& ngrow,
        const int level,
        bool redistribute,
        std::optional<const amrex::Real> initial_value
    )
    {
        // create name
        //     Add the suffix "[level=level]"
        name.append("[level=").append(std::to_string(level)).append("]");

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
        std::string /* other_key */
    )
    {
        throw std::runtime_error("MultiFabRegister::alloc_like not yet implemented");

        // Checks
        // TODO: does the key already exist? error
    }

    /** title
     *
     * body body
     * body
     *
     * @param name ...
     * @return ...
     */
    amrex::MultiFab*
    MultiFabRegister::get (
        std::string name
    )
    {
        if (m_mf_register.count(name) == 0) {
            throw std::runtime_error("MultiFabRegister::get name does not exist in register: " + name);
        }
        amrex::MultiFab & mf = m_mf_register[name].m_mf;

        return &mf;
    }

    /** title
     *
     * body body
     * body
     *
     * @return ...
     */
    std::vector<std::string>
    MultiFabRegister::list ()
    {
        std::vector<std::string> names;
        names.reserve(m_mf_register.size());
        for (auto const & str : m_mf_register) { names.push_back(str.first); }

        return names;
    }


    /** title
     *
     * body body
     * body
     *
     * @param name ...
     * @return ...
     */
    void
    MultiFabRegister::erase (
        std::string name
    )
    {
        if (m_mf_register.count(name) != 1) {
            throw std::runtime_error("MultiFabRegister::remove name does not exist in register: " + name);
        }
        m_mf_register.erase(name);
    }


} // namespace ablastr::fields
