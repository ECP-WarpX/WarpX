/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "HardEdgedQuadrupole.H"
#include "Utils/WarpXUtil.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

HardEdgedQuadrupole::HardEdgedQuadrupole ()
{

    using namespace amrex::literals;

    amrex::ParmParse pp_element_name("quad");

    // All elements should have this block of code that sets up the z location information
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
    d_ze.resize(zcenters.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zcenters.begin(), zcenters.end(), d_zcenters.begin());

    // This is specific to the element type
    amrex::Vector<amrex::Real> dEdx;
    if (queryArrWithParser(pp_element_name, "dEdx", dEdx)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zs.size() == dEdx.size(),
                     "quad: dEdx must have the same length and the zstarts and zends");
    } else {
        dEdx.resize(nelements, 0._rt);
    }

    amrex::Vector<amrex::Real> dBdx;
    if (queryArrWithParser(pp_element_name, "dBdx", dBdx)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zs.size() == dBdx.size(),
                     "quad: dBdx must have the same length and the zstarts and zends");
    } else {
        dBdx.resize(nelements, 0._rt);
    }

    d_dEdx.resize(dEdx.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dEdx.begin(), dEdx.end(), d_dEdx.begin());
    d_dBdx.resize(dBdx.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dBdx.begin(), dBdx.end(), d_dBdx.begin());

}

