/* Copyright 2019-2022 Andrew Myers, Luca Fedeli, Maxence Thevenet,
 *                     Weiqun Zhang, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_CONSTANT_H_
#define ABLASTR_CONSTANT_H_

#include <AMReX_REAL.H>


/** Numerical compile-time constants */
namespace ablastr::constant
{
    /** Mathematical constants */
    namespace math
    {
        //! ratio of a circle's circumference to its diameter
        static constexpr amrex::Real pi = static_cast<amrex::Real>(3.14159265358979323846);

        //! https://tauday.com/tau-manifesto
        static constexpr amrex::Real tau = static_cast<amrex::Real>(2.0 * 3.14159265358979323846);
    } // namespace math

    /** Physical constants
     *
     * Values are the 2018 CODATA recommended values
     * https://physics.nist.gov/cuu/Constants/index.html
     *
     * New additions here should also be considered for addition to
     * `warpx_constants` in WarpXUtil.cpp's `makeParser`, so that they're
     * available in parsing and evaluation of PICMI expressions, as well
     * as the corresponding Python definitions
     */
    namespace SI
    {
        //! vacuum speed of light [m/s]
        static constexpr auto c = static_cast<amrex::Real>( 299'792'458. );
        //! vacuum permittivity: dielectric permittivity of vacuum [F/m]
        static constexpr auto ep0 = static_cast<amrex::Real>( 8.8541878128e-12 );
        //! vacuum permeability: magnetic permeability of vacuum = 4.0e-7 * pi [H/m]
        static constexpr auto mu0 = static_cast<amrex::Real>( 1.25663706212e-06 );
        //! elementary charge [C]
        static constexpr auto q_e = static_cast<amrex::Real>( 1.602176634e-19 );
        //! electron mass [kg]
        static constexpr auto m_e = static_cast<amrex::Real>( 9.1093837015e-31 );
        //! proton mass [kg]
        static constexpr auto m_p = static_cast<amrex::Real>( 1.67262192369e-27 );
        //! dalton: unified atomic mass unit [kg]
        static constexpr auto m_u = static_cast<amrex::Real>( 1.66053906660e-27 );

        //! reduced Planck Constant = h / tau [J*s]
        static constexpr auto hbar = static_cast<amrex::Real>( 1.054571817e-34 );
        //! mu0/(4*MathConst::pi)*q_e*q_e*c/hbar
        static constexpr auto alpha = static_cast<amrex::Real>( 0.007297352573748943 );
        //! 1./(4*MathConst::pi*ep0) * q_e*q_e/(m_e*c*c)
        static constexpr auto r_e = static_cast<amrex::Real>( 2.817940326204929e-15 );
        //! (2.*alpha*alpha*ep0*ep0*hbar*hbar*hbar)/(45.*m_e*m_e*m_e*m_e*c*c*c*c*c)
        static constexpr double xi = 1.3050122447005176e-52;
        //! This should be usable for single precision instead of xi; very close to smallest float32 number possible (1.2e-38)
        static constexpr auto xi_c2 = static_cast<amrex::Real>( 1.1728865132395492e-35 );

        //! Boltzmann constant (exact) [J/K]
        static constexpr auto kb = static_cast<amrex::Real>( 1.380649e-23 );
    } // namespace SI
} // namespace ablastr::constant

#endif // ABLASTR_CONSTANT_H_
