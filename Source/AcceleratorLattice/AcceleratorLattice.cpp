/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "AcceleratorLattice.H"
#include "LatticeElements/Drift.H"
#include "LatticeElements/HardEdgedQuadrupole.H"
#include "LatticeElements/HardEdgedPlasmaLens.H"

#include <AMReX_REAL.H>

AcceleratorLattice::AcceleratorLattice ()
{

    using namespace amrex::literals;

    /* Get the inputs for and initialize all of the lattice element types */
    amrex::ParticleReal z_location = 0._prt;
    ReadLattice("lattice", z_location);

    h_quad.WriteToDevice();
    h_plasmalens.WriteToDevice();
}

void
AcceleratorLattice::ReadLattice (std::string const & root_name, amrex::ParticleReal & z_location)
{
    amrex::ParmParse pp_lattice(root_name);
    std::vector<std::string> lattice_elements;
    pp_lattice.queryarr("elements", lattice_elements);

    if (!lattice_elements.empty()) {
        m_lattice_defined = true;
    }

    // Loop through lattice elements
    for (std::string const & element_name : lattice_elements) {
        // Check the element type
        amrex::ParmParse pp_element(element_name);
        std::string element_type;
        pp_element.get("type", element_type);

        // Initialize the corresponding element according to its type
        if (element_type == "drift") {
            h_drift.AddElement(pp_element, z_location);
        }
        else if (element_type == "quad") {
            h_quad.AddElement(pp_element, z_location);
        }
        else if (element_type == "plasmalens") {
            h_plasmalens.AddElement(pp_element, z_location);
        }
        else if (element_type == "lattice") {
            ReadLattice(element_name, z_location);
        }
        else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                "ERROR: Unknown accelerator lattice element type " + element_type);
        }
    }
}

void
AcceleratorLattice::InitElementFinder (int const lev, amrex::BoxArray const & ba, amrex::DistributionMapping const & dm)
{
    if (m_lattice_defined) {
        m_element_finder = std::make_unique<amrex::LayoutData<LatticeElementFinder>>(ba, dm);
        for (amrex::MFIter mfi(*m_element_finder); mfi.isValid(); ++mfi)
        {
            (*m_element_finder)[mfi].InitElementFinder(lev, mfi, *this);
        }
    }
}

void
AcceleratorLattice::UpdateElementFinder (int const lev)
{
    if (m_lattice_defined) {
        for (amrex::MFIter mfi(*m_element_finder); mfi.isValid(); ++mfi)
        {
            (*m_element_finder)[mfi].UpdateIndices(lev, mfi, *this);
        }
    }
}

LatticeElementFinderDevice
AcceleratorLattice::GetFinderDeviceInstance (WarpXParIter const& a_pti, int const a_offset) const
{
    LatticeElementFinder & finder = (*m_element_finder)[a_pti];
    return finder.GetFinderDeviceInstance(a_pti, a_offset, *this);
}
