/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"
#include "Utils/WarpX_Complex.H"


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
        if (ppc.query(it->c_str(), v)) {
            m_parser.setConstant(*it, v);
            it = symbols.erase(it);
        } else {
            ++it;
        }
    }
    for (auto const& s : symbols) { // make sure there no unknown symbols
        amrex::Abort("Laser Profile: Unknown symbol "+s);
    }

    m_gpu_parser = std::make_unique< ParserWrapper<3> >(m_parser);
}

void
WarpXLaserProfiles::FieldFunctionLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    auto parser = getParser(m_gpu_parser);
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        amplitude[i] = parser(Xp[i], Yp[i], t);
    });
}
