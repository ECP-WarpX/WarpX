/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ELECTROSTATICSOLVER_H_
#define ELECTROSTATICSOLVER_H_

#include "Utils/WarpXUtil.H"

#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

namespace ElectrostaticSolver {

struct PhiCalculatorEB {

    amrex::Real t;
    amrex::ParserExecutor<4> potential_eb;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (const amrex::Real x, const amrex::Real z) const noexcept
    {
        return potential_eb(x, 0.0, z, t);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (const amrex::Real x, const amrex::Real y, const amrex::Real z) const noexcept
    {
        return potential_eb(x, y, z, t);
    }
};

class BoundaryValueHandler {

public:
    std::string potential_eb_str = "0";
    // TODO add handling of domain boundaries as well

    void buildParsers ()
    {
        potential_eb_parser = makeParser(potential_eb_str, {"x", "y", "z", "t"});

        potential_eb = potential_eb_parser.compile<4>();
    }

    PhiCalculatorEB getPhiEB (amrex::Real t) const noexcept
    {
        return PhiCalculatorEB{t, potential_eb};
    }

    amrex::ParserExecutor<4> potential_eb;

private:

    amrex::Parser potential_eb_parser;

};

}
#endif
