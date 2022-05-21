/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "HardEdgedPlasmaLens.H"
#include "Utils/WarpXUtil.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

HardEdgedPlasmaLens::HardEdgedPlasmaLens ()
    : LatticeElementBase("plasmalens")
{

    using namespace amrex::literals;

    if (nelements == 0) return;

    amrex::ParmParse pp_element_name("lattice." + m_element_name);

    amrex::Vector<amrex::Real> dEdx;
    if (queryArrWithParser(pp_element_name, "dEdx", dEdx)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(nelements == dEdx.size(),
                     "plasmalens: dEdx must have the same length and the zstarts and zends");
    } else {
        dEdx.resize(nelements, 0._rt);
    }

    amrex::Vector<amrex::Real> dBdx;
    if (queryArrWithParser(pp_element_name, "dBdx", dBdx)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(nelements == dBdx.size(),
                     "plasmalens: dBdx must have the same length and the zstarts and zends");
    } else {
        dBdx.resize(nelements, 0._rt);
    }

    d_dEdx.resize(dEdx.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dEdx.begin(), dEdx.end(), d_dEdx.begin());
    d_dBdx.resize(dBdx.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dBdx.begin(), dBdx.end(), d_dBdx.begin());

}

HardEdgedPlasmaLensDevice
HardEdgedPlasmaLens::GetDeviceInstance () const
{
    HardEdgedPlasmaLensDevice result;
    result.InitHardEdgedPlasmaLensDevice(*this);
    return result;
}

void
HardEdgedPlasmaLensDevice::InitHardEdgedPlasmaLensDevice (HardEdgedPlasmaLens const& h_plasmalens)
{

    nelements = h_plasmalens.nelements;

    if (nelements == 0) return;

    d_zs_arr = h_plasmalens.d_zs.data();
    d_ze_arr = h_plasmalens.d_ze.data();

    d_dEdx_arr = h_plasmalens.d_dEdx.data();
    d_dBdx_arr = h_plasmalens.d_dBdx.data();

}
