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
    amrex::MultiFab& getsigma_mf  () { return (*m_sigma_mf); }
    /** return MultiFab, epsilon (permittivity) of the medium. */
    amrex::MultiFab& getepsilon_mf  () { return (*m_eps_mf); }
    /** return MultiFab, mu (permeability) of the medium. */
    amrex::MultiFab& getmu_mf  () { return (*m_mu_mf); }
    /** return MultiFab */
    amrex::iMultiFab& get_magnetic_material_id_mf  () { return (*m_magnetic_material_id_mf); }

    /** Bool flag on whether there are any magnetic materials present. */
    bool is_magnetic_material_present () { return m_magnetic_material_present_flag; }

    /** Initializes the Multifabs storing macroscopic properties
     *  with user-defined functions(x,y,z). */
     void InitializeMacroMultiFabUsingParser (amrex::MultiFab *macro_mf,
                                  amrex::ParserExecutor<3> const& macro_parser,
                                  int lev);

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
    /** iMultifab for magnetic objects, where the integer refers specifies to the ID 
     * of which magnetic material is at the indexed location. The properties of each
     * magnetic material are storing in m_magnetic_objects. Note that -1 corresponds to 
     * non-magnetic and defaults to m_mu_mf. */
    std::unique_ptr<amrex::iMultiFab> m_magnetic_material_id_mf;

    /** Stores initialization type for conductivity : constant or parser */
    std::string m_sigma_s = "constant";
    /** Stores initialization type for permittivity : constant or parser */
    std::string m_epsilon_s = "constant";
    /** Stores initialization type for permeability : constant or parser */
    std::string m_mu_s = "constant";

    /** string for storing parser function for conductivity */
    std::string m_str_sigma_function;
    /** string for storing parser function for permittivity */
    std::string m_str_epsilon_function;
    /** string for storing parser function for permeability */
    std::string m_str_mu_function;
    /** string for storing parser function for the integer ID of magnetic material.
     * This index specifies which element to use in m_magnetic_materials. If the ID is -1,
     * this indicates an absence of any magnetic material present at that location. */
    std::string m_str_magnetic_material_id_function;

    /** Parser Wrapper for conductivity */
    std::unique_ptr<amrex::Parser> m_sigma_parser;
    /** Parser Wrapper for permittivity */
    std::unique_ptr<amrex::Parser> m_epsilon_parser;
    /** Parser Wrapper for permeability */
    std::unique_ptr<amrex::Parser> m_mu_parser;
    /** Parser Wrapper for the integer ID of magnetic material.
     * This index specifies which element to use in m_magnetic_materials. If the ID is -1,
     * this indicates an absence of any magnetic material present at that location. */
    std::unique_ptr<amrex::Parser> m_magnetic_material_id_parser;

    /** Vector for storing magnetic materials present */
    amrex::Vector<MagneticMaterial> m_magnetic_materials;

    /** Bool flag on whether magnetic material is present. */
    bool m_magnetic_material_present_flag = false;

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
         const amrex::Real fac1 = 0.5_rt * sigma * dt / epsilon;
         const amrex::Real alpha = (1._rt - fac1)/(1._rt + fac1);
         return alpha;
     }

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real beta (amrex::Real const sigma,
                              amrex::Real const epsilon,
                              amrex::Real dt) {
         using namespace amrex;
         const amrex::Real fac1 = 0.5_rt * sigma * dt / epsilon;
         const amrex::Real beta = dt / ( epsilon * (1._rt + fac1) );
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
        const amrex::Real fac1 = sigma * dt / epsilon;
        const amrex::Real alpha = (1._rt)/(1._rt + fac1);
        return alpha;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real beta (amrex::Real const sigma,
                            amrex::Real const epsilon,
                            amrex::Real dt) {
        using namespace amrex;
        const amrex::Real fac1 = sigma * dt / epsilon;
        const amrex::Real beta = dt / ( epsilon * (1._rt + fac1) );
        return beta;
    }
};

// class MagneticObject {

// public:


// private:
//     /** String for storing parser function, which specifies where the object is location. */
//     std::string m_str_material_location_function;
//     /** Location parser wrappers. Describes where the object is location */
//     std::unique_ptr<amrex::Parser> m_material_location_parser;
//     /** Specifies the hysteretic properties of the magnetic object using the Jiles-Atherton model */
//     JAModelParameters m_ja_model_parameters;
// };

class MagneticMaterial {

public:

    void ComputeH(
        amrex::RealVect H,
        amrex::RealVect B_at_t_plus_dt_over_2,
        amrex::RealVect B_at_t_minus_dt_over_2);

};

struct MagneticMaterials {

     static const JAModelParameters M330_50_A_steel;

};

const JAModelParameters MagneticMaterials::M330_50_A_steel(7.45E-5_rt /* alpha (unitless) */,
                                                           26.1_rt    /* a (A/m) */,
                                                           1.26E6_rt  /* Ms (A/m) */,
                                                           52.3_rt    /* k (A/m) */,
                                                           0.13_rt    /* c (unitless) */);

#endif // WARPX_MACROSCOPIC_PROPERTIES_H_
