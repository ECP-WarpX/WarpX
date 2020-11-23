/* Copyright 2019 Andrew Myers, Luca Fedeli, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_CONST_H_
#define WARPX_CONST_H_

#include <AMReX_REAL.H>
#include <cmath>

// Math constants
namespace MathConst
{
    static constexpr amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);
}

// Physical constants. Values are the 2018 CODATA recommended values
// https://physics.nist.gov/cuu/Constants/index.html
namespace PhysConst
{
    static constexpr auto c     = static_cast<amrex::Real>( 299'792'458. );
    static constexpr auto ep0   = static_cast<amrex::Real>( 8.8541878128e-12 );
    static constexpr auto mu0   = static_cast<amrex::Real>( 1.25663706212e-06 );
    static constexpr auto q_e   = static_cast<amrex::Real>( 1.602176634e-19 );
    static constexpr auto m_e   = static_cast<amrex::Real>( 9.1093837015e-31 );
    static constexpr auto m_p   = static_cast<amrex::Real>( 1.67262192369e-27 );
    static constexpr auto hbar  = static_cast<amrex::Real>( 1.054571817e-34 );
    static constexpr auto alpha = static_cast<amrex::Real>( 0.007297352573748943 );//mu0/(4*MathConst::pi)*q_e*q_e*c/hbar;
    static constexpr auto r_e   = static_cast<amrex::Real>( 2.817940326204929e-15 );//1./(4*MathConst::pi*ep0) * q_e*q_e/(m_e*c*c);
    static constexpr double xi  = 1.3050122447005176e-52; //(2.*alpha*alpha*ep0*ep0*hbar*hbar*hbar)/(45.*m_e*m_e*m_e*m_e*c*c*c*c*c);
    static constexpr auto xi_c2 = static_cast<amrex::Real>( 1.1728865132395492e-35 ); // This should be usable for single precision, though
    // very close to smallest number possible (1.2e-38)
}

#endif
