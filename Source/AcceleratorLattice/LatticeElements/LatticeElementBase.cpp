/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "LatticeElementBase.H"
#include "Utils/WarpXUtil.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

LatticeElementBase::LatticeElementBase (std::string const& element_name)
{
    using namespace amrex::literals;

    m_element_name = element_name;
    amrex::ParmParse pp_element_name("lattice." + m_element_name);

    amrex::Vector<amrex::Real> zs;
    amrex::Vector<amrex::Real> ze;

    queryArrWithParser(pp_element_name, "zstarts", zs);
    queryArrWithParser(pp_element_name, "zends", ze);

    nelements = static_cast<int>(zs.size());

    if (nelements == 0) return;

    CheckElementCorrectness(zs, ze);

    d_zs.resize(zs.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zs.begin(), zs.end(), d_zs.begin());
    d_ze.resize(ze.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ze.begin(), ze.end(), d_ze.begin());

}

void
LatticeElementBase::CheckElementCorrectness (amrex::Vector<amrex::Real> const & zs, amrex::Vector<amrex::Real> const & ze)
{

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(zs.size() == ze.size(),
                 m_element_name + ": zstarts must have the same length and zends");

    // Make sure that elements have ze > zs
    for (int i = 0 ; i < static_cast<int>(zs.size()) ; i++) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(ze[i] > zs[i],
            m_element_name + ": the zends must be greater than the zstarts, but for element "
              + std::to_string(i) + ", ze=" + std::to_string(ze[i]) +
              " is not greater than zs=" + std::to_string(zs[i]) + ".");
    }

    // Make sure that the elements are in increasing order
    for (int i = 1 ; i < static_cast<int>(zs.size()) ; i++) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(zs[i] > ze[i-1],
            m_element_name + ": the elements must not overlap, but for element "
              + std::to_string(i) + ", zs=" + std::to_string(zs[i]) + " is not greater than ze="
              + std::to_string(ze[i-1]) + " of the previous element.");
    }

}
