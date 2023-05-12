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

     /** return MultiFab, sigma (conductivity) of the medium. */
     amrex::MultiFab& getsigma_mf  () {return (*m_sigma_mf);}
     /** return MultiFab, epsilon (permittivity) of the medium. */
     amrex::MultiFab& getepsilon_mf  () {return (*m_eps_mf);}
     /** return MultiFab, mu (permeability) of the medium. */
     amrex::MultiFab& getmu_mf  () {return (*m_mu_mf);}
     /** return MultiFab */
     amrex::iMultiFab& getferromagnetic_id_mf  () {return (*m_ferromagnetic_id_mf);}

     /** Bool flag on whether there are any ferromagnetic materials present. */
     bool is_ferromagnetic_material_present () { return m_ferromagnetic_material_present_flag; }

     /** Initializes the Multifabs storing macroscopic properties
      *  with user-defined functions(x,y,z).
      */
     void InitializeMacroMultiFabUsingParser (amrex::MultiFab *macro_mf,
                                  amrex::ParserExecutor<3> const& macro_parser,
                                  const int lev);

     /** Gpu Vector with index type of the conductivity multifab */
     amrex::GpuArray<int, 3> sigma_IndexType;
     /** Gpu Vector with index type of the permittivity multifab */
     amrex::GpuArray<int, 3> epsilon_IndexType;
     /** Gpu Vector with index type of the permeability multifab */
     amrex::GpuArray<int, 3> mu_IndexType;
     /** Gpu Vector with index type of the Ex multifab */
     amrex::GpuArray<int, 3> Ex_IndexType;
     /** Gpu Vector with index type of the Ey multifab */
     amrex::GpuArray<int, 3> Ey_IndexType;
     /** Gpu Vector with index type of the Ez multifab */
     amrex::GpuArray<int, 3> Ez_IndexType;
     /** Gpu Vector with index type of coarsening ratio with default value (1,1,1) */
     amrex::GpuArray<int, 3> macro_cr_ratio;

private:

     /** Conductivity, sigma, of the medium */
     amrex::Real m_sigma = 0.0;
     /** Permittivity, epsilon, of the medium */
     amrex::Real m_epsilon = PhysConst::ep0;
     /** Permeability, mu, of the medium */
     amrex::Real m_mu = PhysConst::mu0;
     /** Multifab for m_sigma */
     std::unique_ptr<amrex::MultiFab> m_sigma_mf;
     /** Multifab for m_epsilon */
     std::unique_ptr<amrex::MultiFab> m_eps_mf;
     /** Multifab for m_mu */
     std::unique_ptr<amrex::MultiFab> m_mu_mf;
     /** iMultifab for m_ferromagnetic_id. The integer refers specifies which ferromagnetic material is at
      *  the indexed location. Note that 0 (or maybe -1 ...) corresponds to non-ferromagnetic and defaults to m_mu. 
      */
     std::unique_ptr<amrex::iMultiFab> m_ferromagnetic_id_mf;

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

     /** Vector for storing ferromagnetic objects */
     std::vector<FerromagneticObject> m_ferro_objects;

     /** Bool flag on whether ferromagnetic material is present. */
     bool m_ferromagnetic_material_present_flag = false;




};

class FerromagneticObject {

private:
     /** String for storing parser function, which specifies where the object is location. */
     std::string m_str_ferro_location_function;
     /** Location parser wrappers. Describes where the object is location */
     std::unique_ptr<amrex::Parser> m_ferro_location_parser;
     /** Specifies the hysteretic properties of the ferromagnetic object using the Jiles-Atherton model */
     JAModelParameters m_ferro_ja_model_properties;
};

struct FerromagneticMaterials {
public:
     static JAModelParameters M330_50_A_steel JAModelParameters(alphaj = 7.45E-5_rt,
                                                                a = 26.1_rt /* A/m */,
                                                                Ms = 1.26E6_rt /* A/m */,
                                                                k = 52.3_rt /* A/m */,
                                                                c = 0.13_rt);
}

/**
 * \brief This struct contains the five parameters necessary to describe ferromagnetic
 * hysteresis using the Jiles-Atherton (JA) model. 
 */
struct JAModelParameters {
     
public:

     JAModelParameters (amrex::Real alphaj, 
                        amrex::Real a,
                        amrex::Real Ms,
                        amrex::Real k,
                        amrex::Real c) {
          m_alphaj = alphaj;
          m_a = a;
          m_Ms = Ms;
          m_k = k;
          m_c = c;
     }
     /** The `alpha_j` parameter of the JA model is used to define the effective field `H_e = H + alpha_j M`, 
      *  where `H` is the magnetic field and `M` is the total magnetization. It is a unitless parameter. */
     amrex::Real alphaj () { return m_alphaj; }
     /** The `a` parameter of the JA model describes the scale of 
      *  the `H` field at which the ferromagnetic material saturates.
      *  It is used in the Langevin function as `L(H_e / a)`, and
      *  has same units as the H field (A/m).
      */
     amrex::Real a () { return m_a; }
     /** The `M_s` parameter of the JA model specifies the magnetization 
      *  of the ferromagnetic material at maximum saturation.
      *  It has same units as the `H` field (A/m).
      */
     amrex::Real Msaturation () { return m_Ms; }
     /** The `k` parameter of the JA model describes the scale by which 
      *  the irreversible magnetization `M_irr` changes due to a change 
      *  in the effective field `H_e`. In particular, 
      *       `dM_irr / dH_e = (M_an - M_irr) / k`, 
      *  where `M_an` is the anhysteretic magnetization.
      *  It has the same units as the `H` field (A/m). 
      */
     amrex::Real k () { return m_k; }
     /** The `c` parameter of the JA is used in the equation 
      *        `M_rev = c (M_an - M_irr)`, 
      *  where `M_rev`, `M_an`, and `M_irr` are the reversible, anhysteretic, 
      *  and irreversible magnetizations, respectively. 
      *  This equation corresponds to the bulge of the walls between magnetic
      *  domains for small displacement. Given that the total magnetization `M`
      *  is the sum of the reversible and irreversible parts, we can write 
      *  `M = c M_an + (1-c) M_irr`. Note `c` is a unitless parameter. 
      */
     amrex::Real c () { return m_c; }
     
private:
     amrex::Real m_alphaj;
     amrex::Real m_a;
     amrex::Real m_Ms;
     amrex::Real m_k;
     amrex::Real m_c;
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
