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
    amrex::iMultiFab& getferro_material_id_mf  () {return (*m_ferro_material_id_mf);}

    /** Bool flag on whether there are any ferromagnetic materials present. */
    bool is_ferromagnetic_material_present () { return m_ferromagnetic_material_present_flag; }

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
    /** iMultifab for ferromagnetic objects, where the integer refers specifies to the ID 
     * of which ferromagnetic material is at the indexed location. The properties of each
     * ferromagnetic material are storing in m_ferro_objects. Note that -1 corresponds to 
     * non-ferromagnetic and defaults to m_mu_mf. */
    std::unique_ptr<amrex::iMultiFab> m_ferro_material_id_mf;

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
    /** string for storing parser function for the integer ID of ferromagnetic material.
     * This index specifies which element to use in m_ferro_materials. If the ID is -1,
     * this indicates an absence of any ferromagnetic material present at that location. */
    std::string m_str_ferro_material_id_function;

    /** Parser Wrapper for conductivity */
    std::unique_ptr<amrex::Parser> m_sigma_parser;
    /** Parser Wrapper for permittivity */
    std::unique_ptr<amrex::Parser> m_epsilon_parser;
    /** Parser Wrapper for permeability */
    std::unique_ptr<amrex::Parser> m_mu_parser;
    /** Parser Wrapper for the integer ID of ferromagnetic material.
     * This index specifies which element to use in m_ferro_materials. If the ID is -1,
     * this indicates an absence of any ferromagnetic material present at that location. */
    std::unique_ptr<amrex::Parser> m_ferro_material_id_parser;

    /** Vector for storing ferromagnetic materials present */
    amrex::Vector<FerromagneticMaterial> m_ferro_materials;

    /** Bool flag on whether ferromagnetic material is present. */
    bool m_ferromagnetic_material_present_flag = false;

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

<<<<<<< HEAD
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
=======
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
>>>>>>> 7eada29a6981a06651a90a4e46ca0eed1a2e7676

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

<<<<<<< HEAD
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
=======
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
>>>>>>> 7eada29a6981a06651a90a4e46ca0eed1a2e7676

};

class FerromagneticMaterial {

public:

    void ComputeH(
        amrex::RealVect H,
        amrex::RealVect B_at_t_plus_dt_over_2,
        amrex::RealVect B_at_t_minus_dt_over_2);

};

// class FerromagneticObject {

// public:


// private:
//     /** String for storing parser function, which specifies where the object is location. */
//     std::string m_str_ferro_location_function;
//     /** Location parser wrappers. Describes where the object is location */
//     std::unique_ptr<amrex::Parser> m_ferro_location_parser;
//     /** Specifies the hysteretic properties of the ferromagnetic object using the Jiles-Atherton model */
//     JAModelParameters m_ja_model_parameters;
// };

/**
 * \brief This struct contains the five parameters necessary to describe magnetic
 * hysteresis of a material using the Jiles-Atherton (JA) model.
 */
struct JAModelParameters {

public:
    /**
     * \brief Specify the parameters of the Jiles-Atherton (JA) model for a magnetic material.
     * 
     * \param alpha 
     * (unitless) Defines the effective field `H_e = H + alpha M`,
     * where `H` is the magnetic field and `M` is the total magnetization.
     * Quantifies interdomain coupling in the magnetic material.
     * 
     * \param a 
     * (A/m) Describes the scale of the `H` field at which the magnetic 
     * material saturates. Used in the Langevin function as `L(H_e / a)`.
     * Quantifies domain walls density in the magnetic material.
     * 
     * \param Ms
     * (A/m) Saturation magnetization of the magnetic material.
     * 
     * \param k 
     * (A/m) Quantifies average energy required to break a pinning site in the magnetic material.
     * Describes the scale by which the irreversible magnetization `M_irr` changes due to a change
     * in the effective field `H_e`. In particular, `dM_irr / dH_e = (M_an - M_irr) / k`,
     *  where `M_an` is the anhysteretic magnetization.
     * 
     * \param c 
     * (unitless) Quantifies magnetization reversiblity: `M_rev = c (M_an - M_irr)`,
     * where `M_rev`, `M_an`, and `M_irr` are the reversible, anhysteretic, and irreversible
     * magnetizations, respectively. This equation corresponds to the bulge of the domain walls 
     * for small displacement. The total magnetization is
     * `M = M_rev + M_irr = c M_an + (1-c) M_irr`.
     */
    JAModelParameters (amrex::Real a_alpha, amrex::Real a_a, amrex::Real a_Ms, amrex::Real a_k, amrex::Real a_c) {
        m_alpha = a_alpha;
        m_a = a_a;
        m_Ms = a_Ms;
        m_k = a_k;
        m_c = a_c;
        m_Ms_over_a = a_Ms / a_a;
        m_c_times_Ms_over_a = a_c * a_Ms / a_a;
    }
    /**(unitless) Defines the effective field `H_e = H + alpha M`,
     * where `H` is the magnetic field and `M` is the total magnetization.
     * Quantifies interdomain coupling in the magnetic material. */
    amrex::Real alpha () { return m_alpha; }

    /**(A/m) Describes the scale of the `H` field at which the magnetic 
     * material saturates. Used in the Langevin function as `L(H_e / a)`.
     * Quantifies domain walls density in the magnetic material. */
    amrex::Real a () { return m_a; }

    /**(A/m) Saturation magnetization of the magnetic material. */
    amrex::Real Ms () { return m_Ms; }

    /**(A/m) Quantifies average energy required to break a pinning site in the magnetic material.
     * Describes the scale by which the irreversible magnetization `M_irr` changes due to a change
     * in the effective field `H_e`. In particular, `dM_irr / dH_e = (M_an - M_irr) / k`,
     * where `M_an` is the anhysteretic magnetization.
     */
    amrex::Real k () { return m_k; }

    /** (unitless) Quantifies magnetization reversiblity: `M_rev = c (M_an - M_irr)`,
     * where `M_rev`, `M_an`, and `M_irr` are the reversible, anhysteretic, and irreversible
     * magnetizations, respectively. This equation corresponds to the bulge of the domain walls 
     * for small displacement. The total magnetization is 
     * `M = M_rev + M_irr = c M_an + (1-c) M_irr`.
     */
    amrex::Real c () { return m_c; }

    amrex::Real Ms_over_a () { return m_Ms_over_a; }

    amrex::Real c_times_Ms_over_a () { return m_c_times_Ms_over_a; }

    
private:
    amrex::Real m_alpha = 0._rt;
    amrex::Real m_a = 1._rt;
    amrex::Real m_Ms = 1._rt;
    amrex::Real m_k = 1._rt;
    amrex::Real m_c = 0._rt;
    amrex::Real m_one_minus_alpha = 1._rt;
    amrex::Real m_Ms_over_a = 1._rt;
    amrex::Real m_c_times_Ms_over_a = 0._rt;
};

class JAModel {
public:
    /**
     * \brief Update the magnetization for the Jiles-Atherton model
     * 
     * \param[out] M Magnetization vector, will be updated from t-dt/2 to t+dt/2
     * \param[in] B Magnetic flux density B(t-dt/2)
     * \param[in] dB dB = B(t+dt/2) - B(t-dt/2)
     * \param[in] ja_model_parameters Parameters of the Jiles-Atherton model
     */
    void
    UpdateM (amrex::RealVect &M, const amrex::RealVect &B, 
             const amrex::RealVect &dB, JAModelParameters *ja_model_parameters)
    {
        /* 1 - alpha */
        amrex::Real oneMinusAlpha = 1._rt - ja_model_parameters->alpha();
        /* The anhysteretic magnetization */
        amrex::RealVect M_an;
        /* Jacobian matrix chi_an = dM_an/dH_e. The anhysteretic part to chi */
        amrex::Array<amrex::Real, 6> chi_an;
        /* The effective field H_e = H + alpha M = mu0^{-1} B - (1-alpha) M */
        amrex::RealVect H_e = (1._rt / PhysConst::mu0) * B - oneMinusAlpha * M;
        /* Compute M_an and chi_an */ 
        ComputeAnhystereticContribution(M_an, chi_an, H_e, ja_model_parameters);
        /* The change in the effective field */
        amrex::RealVect dH_e;
        amrex::RealVect dB_over_mu0 = (1._rt / PhysConst::mu0) * dB; // mu0^{-1} dB
        amrex::RealVect M_diff = M_an - M; // M_an - M
        amrex::Real M_diff_mag = M_diff.vectorLength(); // |M_an - M|
        amrex::Real denominator = ja_model_parameters->k() * M_diff_mag; // k |M_an - M|
        if (denominator != 0) {
            // Jacobian matrix chi = dM/dH_e = chi_an + chi_irr
            amrex::Array<amrex::Real, 6> chi(chi_an);
            // Add the irreversible contribution chi_irr.
            chi[0] += M_diff[0] * M_diff[0] / denominator; // xx
#if AMREX_SPACEDIM > 1
            chi[1] += M_diff[1] * M_diff[0] / denominator; // xy
            chi[2] += M_diff[1] * M_diff[1] / denominator; // yy
#endif
#if AMREX_SPACEDIM > 2
            chi[3] += M_diff[2] * M_diff[0] / denominator; // xz
            chi[4] += M_diff[2] * M_diff[1] / denominator; // yz
            chi[5] += M_diff[2] * M_diff[2] / denominator; // zz
#endif
            /* Compute dH_e using the total contribution */
            dH_e = Compute_dHe(chi, dB_over_mu0, oneMinusAlpha);
            /* Check the sign of (M_an - M) dot dH_e */
            if (M_diff.dotProduct(dH_e) < 0) {
                /* Use only anhysteretic part of chi if (M_an - M) dot dH_e < 0 */
                dH_e = Compute_dHe(chi_an, dB_over_mu0, oneMinusAlpha);
            }
        }
        else {
            dH_e = Compute_dHe(chi_an, dB_over_mu0, oneMinusAlpha);
        }
        /* Update the magnetization */
        M = M + (1._rt / oneMinusAlpha) * (dB_over_mu0 - dH_e);
    };

private:

    /**
     * \brief Compute the anhysteretic quantities for the JA model.
     * 
     * \param[out] M_an The anhysteretic magnetization
     * \param[out] chi_an The Jacobian matrix dM_an/dH_e.
     *                    Symmetric, with elements ordered as xx, xy, yy, xz, yz, zz.
     * \param[in] H_e The effective magnetic field
     * \param[in] ja_model_parameters Parameters from the Jiles-Atherton model
     */
    void 
    ComputeAnhystereticContribution (amrex::RealVect &M_an,
                                     amrex::Array<amrex::Real, 6> &chi_an,
                                     const amrex::RealVect &H_e,
                                     JAModelParameters* ja_model_parameters) {
        amrex::Real a = ja_model_parameters->a();
        /* c Ms / a */
        amrex::Real c_Ms_over_a = ja_model_parameters->c_times_Ms_over_a();
        amrex::Real H_e_mag = H_e.vectorLength(); // |H_e|
        if (H_e_mag <= 0._rt) { return; }
        amrex::Real x = H_e_mag / a; // x = |H_e| / a
        amrex::Real L_over_x = Langevin::L_over_x(x); // L(x)/x
        amrex::Real dLdx = Langevin::dLdx(x); // L'(x)
        M_an = ((ja_model_parameters->Ms_over_a()) * L_over_x) * H_e;
        amrex::Real delta = dLdx - L_over_x;  // delta = L'(x) - L(x)/x
        amrex::RealVect H_e_hat = H_e / H_e_mag; // Unit vector in direction of H_e
        // chi_an is a symmetric matrix, so only 6 components need to be specified.
        // The order of elements is xx, xy, yy, xz, yz, zz)
        chi_an[0] = c_Ms_over_a * (L_over_x + delta * H_e_hat[0] * H_e_hat[0]);
#if AMREX_SPACEDIM > 1
        chi_an[1] = c_Ms_over_a * (           delta * H_e_hat[1] * H_e_hat[0]);
        chi_an[2] = c_Ms_over_a * (L_over_x + delta * H_e_hat[1] * H_e_hat[1]);
#endif
#if AMREX_SPACEDIM > 2
        chi_an[3] = c_Ms_over_a * (           delta * H_e_hat[2] * H_e_hat[0]);
        chi_an[4] = c_Ms_over_a * (           delta * H_e_hat[2] * H_e_hat[1]);
        chi_an[5] = c_Ms_over_a * (L_over_x + delta * H_e_hat[2] * H_e_hat[2]);
#endif        
    }

    /**
     * \brief Computes change in effective field as dH_e = (dB/dH_e)^{-1} dB.
     * Since B/mu0 = H_e + (1-alpha) M, we have d(B/mu0)/dH_e = I + (1-alpha) chi,
     * where chi = dM/dH_e.
     * 
     * \param chi Components of symmetric matrix chi = dM/dH_e
     * \param dB_over_mu0 mu0^{-1} dB
     * \param oneMinusAlpha 1 - alpha
     * \return amrex::RealVect dH_e = (dB/dH_e)^{-1} dB
     */
    amrex::RealVect
    Compute_dHe (amrex::Array<amrex::Real, 6> chi, RealVect dB_over_mu0, amrex::Real oneMinusAlpha) {
        /* Jacobian matrix mu0^{-1} dB/dHe = I + (1 - alpha) chi */
        Matrix dB_over_mu0_dHe(AMREX_SPACEDIM, amrex::Vector<amrex::Real>(AMREX_SPACEDIM, 0));
        for (int j = 0; j < AMREX_SPACEDIM; ++j){
            int jp1choose2 = j * (j + 1) / 2;
            dB_over_mu0_dHe[j][j] = 1 + oneMinusAlpha * chi[j + jp1choose2];
            for (int i = 0; i < j; ++i){
                dB_over_mu0_dHe[i][j] = oneMinusAlpha * chi[i + jp1choose2];
                dB_over_mu0_dHe[j][i] = dB_over_mu0_dHe[i][j];
            }
        }
        /* Jacobian matrix mu0 dHe/dB */
        Matrix dHe_dB_over_mu0 = MatrixInverse(dB_over_mu0_dHe, AMREX_SPACEDIM);
        /* Could be more efficient by writing a function that solves x = A^{-1} b directly */
        return MatrixVectorProduct(dHe_dB_over_mu0, dB_over_mu0, AMREX_SPACEDIM);
    }

    typedef amrex::Vector<amrex::Vector<amrex::Real>> Matrix;

    Matrix MatrixInverse (Matrix A, int order) {
        // int order = A.size();
        Matrix B(order, amrex::Vector<amrex::Real>(order * 2, 0));

        // copy A into B and create identity matrix
        for (int i = 0; i < order; i++) {
            for (int j = 0; j < order; j++) {
                B[i][j] = A[i][j];
            }
            B[i][i + order] = 1;
        }

        // perform row operations
        for (int i = 0; i < order; i++) {
            // pivot element
            double pivot = B[i][i];
            if (pivot == 0) {
                // matrix is singular
                return Matrix();
            }

            // divide row i by pivot
            for (int j = 0; j < order * 2; j++) {
                B[i][j] /= pivot;
            }

            // eliminate other rows
            for (int j = 0; j < order; j++) {
                if (i != j) {
                    double factor = B[j][i];
                    for (int k = 0; k < order * 2; k++) {
                        B[j][k] -= factor * B[i][k];
                    }
                }
            }
        }

        // extract inverse matrix
        Matrix inverse(order, amrex::Vector<amrex::Real>(order));
        for (int i = 0; i < order; i++) {
            for (int j = 0; j < order; j++) {
                inverse[i][j] = B[i][j + order];
            }
        }
        return inverse;
    }

    Matrix MatrixProduct (Matrix A, Matrix B){
        int n = A.size();
        Matrix product(n, amrex::Vector<amrex::Real>(n, 0));
        for (int i = 0; i < n; ++i){
            for (int j = 0; j < n; ++j){
                for (int k = 0; k < n; ++k){
                    product[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return product;
    }

    amrex::RealVect MatrixVectorProduct (Matrix A, amrex::RealVect v, int order){
        amrex::RealVect result;
        for (int i = 0; i < order; ++i){
            for (int k = 0; k < order; ++k){
                result[i] += A[i][k] * v[k];
            }
        }
        return result;
    }
    
};

struct Langevin {

private:
    // L(x) = (1/3)x + O(x^3)
    static constexpr amrex::Real l1 = 1._rt/3._rt;
    // Slightly different from -1/45 so that L(x) is continuous
    static constexpr amrex::Real l3 = -0.02222201058412696_rt;
    // Slightly different from -1/15 so that dLdx is continuous
    static constexpr amrex::Real dldx2 = -0.06666559177859188_rt;

public:
    /**
     * \brief The Langevin function L(x).
     * 
     * \param x Input (unitless)
     * \return L(x) = coth(x) - 1/x
     */
    static amrex::Real
    L (amrex::Real x) {
        if (x < 0.01_rt) { // Use Taylor expansion for small x
            return l1 * x + l3 * std::pow(x, 3); 
        }
        return 1._rt / std::tanh(x) - 1._rt / x;
     };

    /**
     * \brief The Langevin function L(x) over x.
     * 
     * \param x Input (unitless)
     * \return L(x)/x = coth(x)/x - 1/x^2
     */
     static amrex::Real
     L_over_x (amrex::Real x) {
          if (x < 0.01_rt) { // Use Taylor expansion for small x
               return l1 + l3 * x * x; 
          }
          amrex::Real xinv = 1._rt / x;
          return xinv * (1._rt / std::tanh(x) - xinv );
     };

    /**
     * \brief The derivative of the Langevin function.
     * 
     * \param x Input (unitless)
     * \return L'(x) = - 1/sinh^2(x) + 1/x^2
     */
     static amrex::Real 
     dLdx (amrex::Real x) {
        if (x < 0.01_rt) { // Use Taylor expansion for small x
               return l1 + dldx2 * x * x; 
          }
          amrex::Real y = std::sinh(x);
          return -1._rt / (y * y) + 1._rt / (x * x);
     };
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
