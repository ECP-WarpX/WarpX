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
    m_element_name = element_name;
}

void
LatticeElementBase::AddElementBase (amrex::ParmParse & pp_element, amrex::ParticleReal & z_location)
{
    // Read in the length of the element and save the start and end, and update z_location
    amrex::ParticleReal ds;
    pp_element.get("ds", ds);

    h_zs.push_back(z_location);
    z_location += ds;
    h_ze.push_back(z_location);

    nelements += 1;
}

void
LatticeElementBase::WriteToDeviceBase ()
{
    d_zs.resize(h_zs.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_zs.begin(), h_zs.end(), d_zs.begin());
    d_ze.resize(h_ze.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_ze.begin(), h_ze.end(), d_ze.begin());
}
