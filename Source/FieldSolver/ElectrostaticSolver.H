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
#include <AMReX_MLNodeTensorLaplacian.H>
#if defined(AMREX_USE_EB) || defined(WARPX_DIM_RZ)
#    include <AMReX_MLEBNodeFDLaplacian.H>
#endif
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

class BoundaryHandler {

public:

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
    bool bcs_set = false;
    std::array<bool,AMREX_SPACEDIM*2> dirichlet_flag;
    bool has_non_periodic = false;
    bool phi_EB_only_t = true;

    void definePhiBCs ();
    void buildParsers ();

    PhiCalculatorEB getPhiEB (amrex::Real t) const noexcept
    {
        return PhiCalculatorEB{t, potential_eb};
    }

    // set default potentials to zero in order for current tests to pass
    // but forcing the user to specify a potential might be better
    std::string potential_xlo_str = "0";
    std::string potential_xhi_str = "0";
    std::string potential_ylo_str = "0";
    std::string potential_yhi_str = "0";
    std::string potential_zlo_str = "0";
    std::string potential_zhi_str = "0";
    std::string potential_eb_str = "0";

    amrex::ParserExecutor<1> potential_xlo;
    amrex::ParserExecutor<1> potential_xhi;
    amrex::ParserExecutor<1> potential_ylo;
    amrex::ParserExecutor<1> potential_yhi;
    amrex::ParserExecutor<1> potential_zlo;
    amrex::ParserExecutor<1> potential_zhi;
    amrex::ParserExecutor<1> potential_eb_t;
    amrex::ParserExecutor<4> potential_eb;

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
