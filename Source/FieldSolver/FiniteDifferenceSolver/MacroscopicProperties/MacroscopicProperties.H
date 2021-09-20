/* Copyright 2020 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_MACROSCOPICPROPERTIES_H_
#define WARPX_MACROSCOPICPROPERTIES_H_

#include "MacroscopicProperties_fwd.H"

#include "Utils/WarpXConst.H"

#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>


enum MacroparameterInitType { ConstantValue, ParserFunction};

/**
 * \brief Functor to return macropameter, either constant value, m_value, or
 *        spatially varying scalar value computed using the parser function, m_parser.
 */

struct GetMacroparameter
{
    /* Type of initialization for macroparameter, constant or parser function */
    MacroparameterInitType m_type;
    /* Constant value of the macroparameter. */
    amrex::Real m_value;
    /* Parser funtion of the spatially varying macroparameter*/
    amrex::ParserExecutor<3> m_parser;
    /**
     * \brief Functor call. This method returns the value of the macroparameter,
     *        or property of the medium needed for the macroscopic Maxwell solver,
     *        at a given location (x,y,z) in the domain.
     *
     * @param[in] x x-coordinate of a given location
     * @param[in] y y-coordinate of a given location
     * @param[in] z z-coordinate of a given location
     * @return value of the macroparameter at (x,y,z).
     *         m_value if init-type is constant
     *         m_parser(x,y,z) if init-type is parser function
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator () (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        using namespace amrex::literals;
        if (m_type == ConstantValue)
        {
            return m_value;
        }
        else if (m_type == ParserFunction)
        {
            return m_parser(x,y,z);
        }
        else
        {
            amrex::Abort("macroparameter init type not valid.");
            return 0.;
        }
        return 0.;
    }
};

/**
 * \brief Functor for conductivity, sigma, of the medium.
 */
struct GetSigmaMacroparameter : GetMacroparameter
{
    /** Constructor to store the type of intialization, m_type, and the value or parser function. */
    GetSigmaMacroparameter () noexcept;
};

/**
 * \brief Functor for permeability, mu, of the medium.
 */
struct GetMuMacroparameter : GetMacroparameter
{
    /** Constructor to store the type of intialization, m_type, and the value or parser function. */
    GetMuMacroparameter () noexcept;
};

/**
 * \brief Functor for permittivity, epsilon, of the medium.
 */
struct GetEpsilonMacroparameter : GetMacroparameter
{
    /** Constructor to store the type of intialization, m_type, and the value or parser function. */
    GetEpsilonMacroparameter () noexcept;
};
/**
 * \brief This class contains the macroscopic properties of the medium needed to
 * evaluate macroscopic Maxwell equation.
 */
class
MacroscopicProperties
{
public:
     MacroscopicProperties (); // constructor
     /** Read user-defined macroscopic properties. Called in constructor. */
     void ReadParameters ();
     /** Initialize multifabs storing macroscopic multifabs */
     void InitData ();

     /** Gpu Vector with index type of the Ex multifab */
     amrex::GpuArray<int, 3> Ex_IndexType;
     /** Gpu Vector with index type of the Ey multifab */
     amrex::GpuArray<int, 3> Ey_IndexType;
     /** Gpu Vector with index type of the Ez multifab */
     amrex::GpuArray<int, 3> Ez_IndexType;
     /** Gpu Vector with index type of the Bx multifab */
     amrex::GpuArray<int, 3> Bx_IndexType;
     /** Gpu Vector with index type of the By multifab */
     amrex::GpuArray<int, 3> By_IndexType;
     /** Gpu Vector with index type of the Bz multifab */
     amrex::GpuArray<int, 3> Bz_IndexType;

     /** Conductivity, sigma, of the medium */
     amrex::Real m_sigma = 0.0;
     /** Permittivity, epsilon, of the medium */
     amrex::Real m_epsilon = PhysConst::ep0;
     /** Permeability, mu, of the medium */
     amrex::Real m_mu = PhysConst::mu0;
     /** Stores initialization type for conductivity : constant or parser */
     std::string m_sigma_s = "constant";
     /** Stores initialization type for permittivity : constant or parser */
     std::string m_epsilon_s = "constant";
     /** Stores initialization type for permeability : constant or parser */
     std::string m_mu_s = "constant";

     /** string for storing parser function */
     std::string m_str_sigma_function;
     std::string m_str_epsilon_function;
     std::string m_str_mu_function;
     /** Parser Wrappers */
     std::unique_ptr<amrex::Parser> m_sigma_parser;
     std::unique_ptr<amrex::Parser> m_epsilon_parser;
     std::unique_ptr<amrex::Parser> m_mu_parser;

};

/**
 * \brief
 * This struct contains only static functions to compute the co-efficients for the
 * Lax-Wendroff scheme of macroscopic Maxwells equations using
 * macroscopic properties, namely, conductivity (sigma), permittivity (epsilon).
 * Permeability of the material, mu, is used as (beta/mu) for the E-update
 * defined in MacroscopicEvolveECartesian().
 */
struct LaxWendroffAlgo {

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real alpha (amrex::Real const sigma,
                               amrex::Real const epsilon,
                               amrex::Real dt) {
         using namespace amrex;
         amrex::Real fac1 = 0.5_rt * sigma * dt / epsilon;
         amrex::Real alpha = (1._rt - fac1)/(1._rt + fac1);
         return alpha;
     }

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real beta (amrex::Real const sigma,
                              amrex::Real const epsilon,
                              amrex::Real dt) {
         using namespace amrex;
         amrex::Real fac1 = 0.5_rt * sigma * dt / epsilon;
         amrex::Real beta = dt / ( epsilon * (1._rt + fac1) );
         return beta;
     }

};

/**
 * \brief
 * This struct contains only static functions to compute the co-efficients for the
 * BackwardEuler scheme of macroscopic Maxwells equations using
 * macroscopic properties, namely, conductivity (sigma) and permittivity (epsilon).
 * Permeability of the material, mu, is used as (beta/mu) for the E-update
 * defined in MacroscopicEvolveECartesian().
 */
struct BackwardEulerAlgo {

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real alpha (amrex::Real const sigma,
                               amrex::Real const epsilon,
                               amrex::Real dt) {
         using namespace amrex;
         amrex::Real fac1 = sigma * dt / epsilon;
         amrex::Real alpha = (1._rt)/(1._rt + fac1);
         return alpha;
     }

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real beta (amrex::Real const sigma,
                              amrex::Real const epsilon,
                              amrex::Real dt) {
         using namespace amrex;
         amrex::Real fac1 = sigma * dt / epsilon;
         amrex::Real beta = dt / ( epsilon * (1._rt + fac1) );
         return beta;
     }

};

#endif // WARPX_MACROSCOPIC_PROPERTIES_H_
