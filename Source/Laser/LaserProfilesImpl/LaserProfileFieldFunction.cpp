/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"

#include "Utils/TextMsg.H"
#include "Utils/WarpX_Complex.H"
#include "Utils/WarpXUtil.H"

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <memory>
#include <set>
#include <string>

using namespace amrex;

void
WarpXLaserProfiles::FieldFunctionLaserProfile::init (
    const amrex::ParmParse& ppl,
    const amrex::ParmParse& ppc,
    CommonLaserParameters /*params*/)
{
    // Parse the properties of the parse_field_function profile
    ppl.get("field_function(X,Y,t)", m_params.field_function);
    m_parser.define(m_params.field_function);
    m_parser.registerVariables({"X","Y","t"});

    std::set<std::string> symbols = m_parser.symbols();
    symbols.erase("X");
    symbols.erase("Y");
    symbols.erase("t"); // after removing variables, we are left with constants
    for (auto it = symbols.begin(); it != symbols.end(); ) {
        Real v;
        if (queryWithParser(ppc, it->c_str(), v)) {
            m_parser.setConstant(*it, v);
            it = symbols.erase(it);
        } else {
            ++it;
        }
    }

    std::stringstream ss;
    for (auto const& s : symbols) ss << " " << s;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(symbols.empty(),
        "Laser Profile: Unknown symbols " + ss.str());

}

void
WarpXLaserProfiles::FieldFunctionLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    auto parser = m_parser.compile<3>();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        amplitude[i] = parser(Xp[i], Yp[i], t);
    });
}
