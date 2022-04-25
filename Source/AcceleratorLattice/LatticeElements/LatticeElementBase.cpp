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
    amrex::ParmParse pp_element_name(m_element_name);

    amrex::Vector<amrex::Real> zs;
    amrex::Vector<amrex::Real> ze;

    queryArrWithParser(pp_element_name, "zstarts", zs);
    queryArrWithParser(pp_element_name, "zends", ze);

    nelements = static_cast<int>(zs.size());

    if (nelements == 0) return;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zs.size() == ze.size(),
                 "quad: zstarts must have the same length and zends");

    amrex::Vector<amrex::Real> zcenters;
    zcenters.resize(nelements + 1, std::numeric_limits<amrex::Real>::lowest());
    for (int i = 1 ; i < nelements ; i++) {
        zcenters[i] = 0.5_rt*(ze[i-1] + zs[i]);
    }
    zcenters[nelements] = std::numeric_limits<amrex::Real>::max();

    d_zs.resize(zs.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zs.begin(), zs.end(), d_zs.begin());
    d_ze.resize(ze.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ze.begin(), ze.end(), d_ze.begin());
    d_zcenters.resize(zcenters.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zcenters.begin(), zcenters.end(), d_zcenters.begin());

}
