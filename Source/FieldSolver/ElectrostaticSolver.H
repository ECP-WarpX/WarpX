/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ELECTROSTATICSOLVER_H_
#define ELECTROSTATICSOLVER_H_

#include "Utils/WarpXUtil.H"

#include <AMReX_Array.H>
#ifdef WARPX_DIM_RZ
#    include <AMReX_MLNodeLaplacian.H>
#else
#    include <AMReX_MLNodeTensorLaplacian.H>
#    ifdef AMREX_USE_EB
#        include <AMReX_MLEBNodeFDLaplacian.H>
#    endif
#endif
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

namespace ElectrostaticSolver {

struct PhiCalculatorEB {

    amrex::Real t;
    amrex::Parser potential_eb_parser;
    amrex::ParserExecutor<4> potential_eb = potential_eb_parser.compile<4>();

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

class BoundaryHandler {

public:

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
    std::array<bool,AMREX_SPACEDIM> dirichlet_flag;
    bool has_Dirichlet = false;

    // set default potentials to zero in order for current tests to pass
    // but forcing the user to specify a potential might be better
    std::string potential_xlo_str = "0";
    std::string potential_xhi_str = "0";
    std::string potential_ylo_str = "0";
    std::string potential_yhi_str = "0";
    std::string potential_zlo_str = "0";
    std::string potential_zhi_str = "0";
    std::string potential_eb_str = "0";

    void definePhiBCs ();

    void buildParsers ()
    {
        potential_xlo_parser = makeParser(potential_xlo_str, {"t"});
        potential_xhi_parser = makeParser(potential_xhi_str, {"t"});
        potential_ylo_parser = makeParser(potential_ylo_str, {"t"});
        potential_yhi_parser = makeParser(potential_yhi_str, {"t"});
        potential_zlo_parser = makeParser(potential_zlo_str, {"t"});
        potential_zhi_parser = makeParser(potential_zhi_str, {"t"});
        potential_eb_parser = makeParser(potential_eb_str, {"x", "y", "z", "t"});

        potential_xlo = potential_xlo_parser.compile<1>();
        potential_xhi = potential_xhi_parser.compile<1>();
        potential_ylo = potential_ylo_parser.compile<1>();
        potential_yhi = potential_yhi_parser.compile<1>();
        potential_zlo = potential_zlo_parser.compile<1>();
        potential_zhi = potential_zhi_parser.compile<1>();
    };

    PhiCalculatorEB getPhiEB (amrex::Real t) const noexcept
    {
        return PhiCalculatorEB{t, potential_eb_parser};
    }

    amrex::ParserExecutor<1> potential_xlo;
    amrex::ParserExecutor<1> potential_xhi;
    amrex::ParserExecutor<1> potential_ylo;
    amrex::ParserExecutor<1> potential_yhi;
    amrex::ParserExecutor<1> potential_zlo;
    amrex::ParserExecutor<1> potential_zhi;

private:

    amrex::Parser potential_xlo_parser;
    amrex::Parser potential_xhi_parser;
    amrex::Parser potential_ylo_parser;
    amrex::Parser potential_yhi_parser;
    amrex::Parser potential_zlo_parser;
    amrex::Parser potential_zhi_parser;
    amrex::Parser potential_eb_parser;
};
}
#endif
