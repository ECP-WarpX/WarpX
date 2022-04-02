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
    : LatticeElement(element_name)
{

    amrex::ParmParse pp_element_name(element_name);

    queryArrWithParser(pp_element_name, "dEdx", m_dEdx);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_zstarts.size() == m_dEdx.size(),
                 element_name + ": The dEdx must have the same length and the zstarts and zends");

}

