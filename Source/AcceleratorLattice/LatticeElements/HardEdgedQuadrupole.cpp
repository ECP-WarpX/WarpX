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
    : LatticeElementBase("quad")
{

    using namespace amrex::literals;

    if (nelements == 0) return;

    amrex::ParmParse pp_element_name(m_element_name);

    amrex::Vector<amrex::Real> dEdx;
    if (queryArrWithParser(pp_element_name, "dEdx", dEdx)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nelements == dEdx.size(),
                     "quad: dEdx must have the same length and the zstarts and zends");
    } else {
        dEdx.resize(nelements, 0._rt);
    }

    amrex::Vector<amrex::Real> dBdx;
    if (queryArrWithParser(pp_element_name, "dBdx", dBdx)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nelements == dBdx.size(),
                     "quad: dBdx must have the same length and the zstarts and zends");
    } else {
        dBdx.resize(nelements, 0._rt);
    }

    d_dEdx.resize(dEdx.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dEdx.begin(), dEdx.end(), d_dEdx.begin());
    d_dBdx.resize(dBdx.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dBdx.begin(), dBdx.end(), d_dBdx.begin());

}

HardEdgedQuadrupoleDevice
HardEdgedQuadrupole::GetDeviceInstance () const
{
    HardEdgedQuadrupoleDevice result;
    result.InitHardEdgedQuadrupoleDevice(*this);
    return result;
}

void
HardEdgedQuadrupoleDevice::InitHardEdgedQuadrupoleDevice (HardEdgedQuadrupole const& h_quad)
{

    nelements = h_quad.nelements;

    if (nelements == 0) return;

    d_zs_arr = h_quad.d_zs.data();
    d_ze_arr = h_quad.d_ze.data();

    d_dEdx_arr = h_quad.d_dEdx.data();
    d_dBdx_arr = h_quad.d_dBdx.data();

}
