/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "LatticeElement.H"
#include "Utils/WarpXUtil.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

LatticeElement::LatticeElement (std::string const & element_name)
{

    amrex::ParmParse pp_element_name(element_name);

    queryArrWithParser(pp_element_name, "zstarts", m_zstarts);
    queryArrWithParser(pp_element_name, "zends", m_zends);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_zstarts.size() == m_zends.size(),
                 element_name + ": The zstarts and zends must have the same length");

    nelements = static_cast<int>(m_zstarts.size());

}

