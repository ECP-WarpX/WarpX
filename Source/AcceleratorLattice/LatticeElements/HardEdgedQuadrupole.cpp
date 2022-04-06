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

    amrex::ParmParse pp_element_name("quad");

    amrex::Vector<amrex::Real> zs;
    amrex::Vector<amrex::Real> ze;
    amrex::Vector<amrex::Real> dEdx;

    queryArrWithParser(pp_element_name, "zstarts", zs);
    queryArrWithParser(pp_element_name, "zends", ze);
    queryArrWithParser(pp_element_name, "dEdx", dEdx);

    nelements = static_cast<int>(zs.size());

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zs.size() == ze.size(),
                 "quad: zstarts must have the same length and zends");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zs.size() == dEdx.size(),
                 "quad: dEdx must have the same length and the zstarts and zends");

    d_zs.resize(zs.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zs.begin(), zs.end(), d_zs.begin());
    d_ze.resize(ze.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ze.begin(), ze.end(), d_ze.begin());
    d_dEdx.resize(dEdx.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dEdx.begin(), dEdx.end(), d_dEdx.begin());

}

