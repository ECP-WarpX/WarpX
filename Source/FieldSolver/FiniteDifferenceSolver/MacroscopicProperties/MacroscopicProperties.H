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
     *  with user-defined functions(x,y,z). */
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

class FerromagneticObject {

public:

private:
    /** String for storing parser function, which specifies where the object is location. */
    std::string m_str_ferro_location_function;
    /** Location parser wrappers. Describes where the object is location */
    std::unique_ptr<amrex::Parser> m_ferro_location_parser;
    /** Specifies the hysteretic properties of the ferromagnetic object using the Jiles-Atherton model */
    JAModelParameters m_ja_model_parameters;
};

struct MagneticMaterials {

     static const JAModelParameters M330_50_A_steel;

};

const JAModelParameters MagneticMaterials::M330_50_A_steel(7.45E-5_rt /* alpha (unitless) */,
                                                           26.1_rt    /* a (A/m) */,
                                                           1.26E6_rt  /* Ms (A/m) */,
                                                           52.3_rt    /* k (A/m) */,
                                                           0.13_rt    /* c (unitless) */);



/**
 * \brief This class contains the five parameters necessary to describe magnetic
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
    }
    /**(unitless) Defines the effective field `H_e = H + alpha M`,
     * where `H` is the magnetic field and `M` is the total magnetization.
     * Quantifies interdomain coupling in the magnetic material. */
    amrex::Real alpha () { return m_alpha; }

    /**Describes the scale of the `H` field at which the magnetic 
     * material saturates. Used in the Langevin function as `L(H_e / a)`.
     * Quantifies domain walls density in the magnetic material. */
    amrex::Real a () { return m_a; }

    /** (A/m) Saturation magnetization of the magnetic material. */
    amrex::Real Ms () { return m_Ms; }

    /**(A/m) Quantifies average energy required to break a pinning site in the magnetic material.
     * Describes the scale by which the irreversible magnetization `M_irr` changes due to a change
     * in the effective field `H_e`. In particular, `dM_irr / dH_e = (M_an - M_irr) / k`,
     * where `M_an` is the anhysteretic magnetization.
     */
    amrex::Real k () { return m_k; }

    /**(unitless) Quantifies magnetization reversiblity: `M_rev = c (M_an - M_irr)`,
     * where `M_rev`, `M_an`, and `M_irr` are the reversible, anhysteretic, and irreversible
     * magnetizations, respectively. This equation corresponds to the bulge of the domain walls 
     * for small displacement. The total magnetization is 
     * `M = M_rev + M_irr = c M_an + (1-c) M_irr`.
     */
    amrex::Real c () { return m_c; }
    
private:
    amrex::Real m_alpha = 0_rt;
    amrex::Real m_a = 1_rt;
    amrex::Real m_Ms = 1_rt;
    amrex::Real m_k = 1_rt;
    amrex::Real m_c = 0_rt;
};

class JAModel {

public:

    void ComputeH(amrex::RealVect & H,
                amrex::RealVect const& B,
                amrex::RealVect const& E,
                amrex::Real dt,
                JAModelParameters* ja_model_parameters)
    {
        amrex::Real alpha = ja_model_parameters->alpha();
        amrex::Real k =  ja_model_parameters->k();
        amrex::RealVect M = (B / PhysConst::mu0) - H;
        amrex::RealVect H_e = H + (alpha * M);
        amrex::RealVect M_an; // The anhysteretic magnetization
        std::array<amrex::Real, 6> chi_an; // The anhysteretic contribution to chi

        amrex::Vector<amrex::Real> v1(3);
        amrex::Vector<amrex::Real> v2(3);
        amrex::Vector<amrex::Real> v3 = v1 - v2;

        ComputeAnhystereticContribution(M_an, chi_an, H_e, ja_model_parameters);

        amrex::RealVect M_diff = M_an - M;
        amrex::Real M_diff_mag = M_diff.vectorLength();
        amrex::Real denominator = k * M_diff_mag;

        std::array<amrex::Real, 6> chi; // The total contribution

        // Now account for the irreversible contribution
        // The order of elements is xx, xy, yy, xz, yz, zz)
        chi[0] = chi_an[0] + M_diff[0] * M_diff[0] / denominator;
#if AMREX_SPACEDIM > 1
        chi[1] = chi_an[1] + M_diff[1] * M_diff[0] / denominator;
        chi[2] = chi_an[2] + M_diff[1] * M_diff[1] / denominator;
#endif
#if AMREX_SPACEDIM > 2
        chi[3] = chi_an[3] + M_diff[2] * M_diff[0] / denominator;
        chi[4] = chi_an[4] + M_diff[2] * M_diff[1] / denominator;
        chi[5] = chi_an[5] + M_diff[2] * M_diff[2] / denominator;
#endif

        Matrix chi_matrix(AMREX_SPACEDIM, amrex::Vector<amrex::Real>(AMREX_SPACEDIM, 0));
        Matrix chi_an_matrix(AMREX_SPACEDIM, amrex::Vector<amrex::Real>(AMREX_SPACEDIM, 0));
        for (int j = 0; j < AMREX_SPACEDIM; ++j){
            chi_matrix[j][j] = chi[i * (3 + i) / 2];
            chi_an_matrix[j][j] = chi_an[i * (3 + i) / 2];
            for (int i = 0; i < j; ++i){
                chi_matrix[i][j] = chi[i + j * (j + 1) / 2];
                chi_matrix[j][i] = chi[i + j * (j + 1) / 2];
                chi_an_matrix[i][j] = chi_an[i + j * (j + 1) / 2];
                chi_an_matrix[j][i] = chi_an[i + j * (j + 1) / 2];
            }
        }

        Matrix dHdB = Compute_dHdB(chi_matrix, alpha);

        


        /* Compute dH = (dH/dB)dB = (dH/dB)(-curl E)dt/2 */

        /* Compute dH_e = (1 - alpha)dH + alpha/mu0 dB 
                        = ((1 - alpha)(dH/dB) + alpha/mu0 identity) dB
                        = ((1 - alpha)(dH/dB) + alpha/mu0 identity) (-curl E) dt/2
                        */

        /** Determine sign of (M_an - M) \cdot dH_e.
         *  If negative, recalculate H by replacing chi with chi_an */

    };


private:

    typedef amrex::Vector<amrex::Vector<amrex::Real>> Matrix;

    // typedef std::array<amrex::Real, AMREX_SPACEDIM * AMREX_SPACEDIM> Matrix;

    /**
     * \brief Compute the derivative matrix dH/dB. 
     * 
     * \param chi Matrix elements of dM/dH_e. Always symmetric.
     * \param alpha Parameter from Jiles-Atherton model.
     * \return Matrix with elements dH_j/dB_i, where i (j) is the row (column) index.
     */
    Matrix Compute_dHdB(Matrix chi, amrex::Real alpha){

        /* Compute I - alpha chi */
        int n = chi.size();
        Matrix I_alphachi(n, amrex::Vector<amrex::Real>(n, 0));
        for (int i = 0; i < n; ++i){
            for (int j = 0; j < n; ++j){
                I_alphachi[i][j] = ((i == j) ? 1._rt : 0._rt) - alpha * chi[i][j];
            }
        }

        /* Compute (I - alpha chi)^(-1) */
        Matrix I_alphachi_inv = MatrixInverse(I_alphachi);

        /* Compute dM_j/dH_i = ([I - alpha chi]^(-1) chi)_{ij} */
        Matrix dMdH = MatrixProduct(I_alphachi_inv, chi);

        Matrix dBdH;
        /* Compute dB_j/dH_i = m0 (I + dM/dH) */
        for (int i = 0; i < n; ++i){
            for (int j = 0; j < n; ++j){
                dBdH[i][j] = PhysConst::mu0 * ( (i == j ? 1._rt : 0._rt) + dMdH[i][j]);
            }
        }

        /* Compute dH/dB = [identity + [identity - alpha chi]^(-1) chi ]^(-1)  */
        return MatrixInverse(dBdH);
    }


    /**
     * \brief 
     * 
     * \param[out] M_an The anhysteretic magnetization
     * \param[out] chi_an The matrix of derivatives of M_an w.r.t. H_e. Turns out to be symmetric in this case, so only 6 components are necessary. The elements are ordered xx, xy, yy, xz, yz, zz.
     * \param[in] H_e The effective magnetic field
     * \param[in] ja_model_parameters Parameters from the Jiles-Atherton model
     */
    void ComputeAnhystereticContribution(amrex::RealVect & M_an,
                                        std::array<amrex::Real, 6>& chi_an,
                                        amrex::RealVect const& H_e,
                                        JAModelParameters* ja_model_parameters){
        amrex::Real a = ja_model_parameters->a();
        amrex::Real Ms = ja_model_parameters->Ms();
        amrex::Real c = ja_model_parameters->c();

        amrex::Real c_Ms_over_a = c * Ms / a;
        amrex::Real H_e_mag = H_e.vectorLength();
        if (H_e_mag <= 0._rt){
            return;
        }
        amrex::Real x = H_e_mag / a;
        amrex::Real L_over_x = Langevin::L_over_x(x);
        amrex::Real dLdx = Langevin::dLdx(x);
        
        
        M_an = ((Ms / a) * L_over_x) * H_e;

        amrex::Real delta = dLdx - L_over_x;
        amrex::RealVect H_e_hat = H_e / H_e_mag; // Division errors could be an issue here

        // for (int k = 0; k < AMREX_SPACEDIM; ++k){
        //     int kchoose2 = k * (k + 1) / 2;
        //     for (int i = 0; i <= k; ++i){
        //         chi_an[i + kchoose2] = c_Ms_over_a * (delta * H_e_hat[i] * H_e_hat[k] + L_over_x);
        //     }
        // }

        // chi_an is a symmetric matrix, so only 6 components need to be specified.
        // The order of elements is xx, xy, yy, xz, yz, zz)
        chi_an[0] = c_Ms_over_a * (delta * H_e_hat[0] * H_e_hat[0] + L_over_x);
#if AMREX_SPACEDIM > 1
        chi_an[1] = c_Ms_over_a * (delta * H_e_hat[1] * H_e_hat[0]           );
        chi_an[2] = c_Ms_over_a * (delta * H_e_hat[1] * H_e_hat[1] + L_over_x);
#endif
#if AMREX_SPACEDIM > 2
        chi_an[3] = c_Ms_over_a * (delta * H_e_hat[2] * H_e_hat[0]           );
        chi_an[4] = c_Ms_over_a * (delta * H_e_hat[2] * H_e_hat[1]           );
        chi_an[5] = c_Ms_over_a * (delta * H_e_hat[2] * H_e_hat[2] + L_over_x);
#endif        
    }


    // void SymmetricMatrixInverse(std::array<amrex::Real, 6>& invS, std::array<amrex::Real, 6> const& S){
    //     amrex::Real det = S[]
    // }

    // template <int n>
    // void MatrixInverse(std::array<amrex::Real>& invM, std::array<amrex::Real, 6> const& S){
    //     amrex::Real det = S[]
    // }

    

    Matrix MatrixInverse(Matrix A) {
        int n = A.size();
        Matrix B(n, amrex::Vector<amrex::Real>(n * 2, 0));

        // copy A into B and create identity matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                B[i][j] = A[i][j];
            }
            B[i][i + n] = 1;
        }

        // perform row operations
        for (int i = 0; i < n; i++) {
            // pivot element
            double pivot = B[i][i];
            if (pivot == 0) {
                // matrix is singular
                return Matrix();
            }

            // divide row i by pivot
            for (int j = 0; j < n * 2; j++) {
                B[i][j] /= pivot;
            }

            // eliminate other rows
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double factor = B[j][i];
                    for (int k = 0; k < n * 2; k++) {
                        B[j][k] -= factor * B[i][k];
                    }
                }
            }
        }

        // extract inverse matrix
        Matrix inverse(n, amrex::Vector<amrex::Real>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inverse[i][j] = B[i][j + n];
            }
        }
        return inverse;
    }

    Matrix MatrixProduct(Matrix A, Matrix B){
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

    amrex::RealVect MatrixVectorProduct(Matrix A, amrex::RealVect v, int order){
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

     static amrex::Real const l1 = 1._rt/3._rt; // L(x) = x/3 + O(x^3)
     static amrex::Real const l3 = -0.02222226972489055_rt; // Slightly different from -1/45 so that L(x) is continuous
     static amrex::Real const dldx2 = -0.06666510671493597_rt; // Slightly different from -1/15 so that dLdx is continuous

public:

    /**
     * \brief The Langevin function L(x).
     * 
     * \param x Input (unitless)
     * \return L(x) = coth(x) - 1/x
     */
     static amrex::Real L(amrex::Real x){
          if (x < 0.01_rt){
               return l1 * x + l3 * pow(x, 3); // Use Taylor expansion for small x
          }
          amrex::Real y =  exp(-2._rt * x);
          return (1._rt + y) / (1._rt - y) - 1._rt / x;
     };

    /**
     * \brief The Langevin function L(x) over x.
     * 
     * \param x Input (unitless)
     * \return L(x)/x = coth(x)/x - 1/x^2
     */
     static amrex::Real L_over_x(amrex::Real x){
          if (x < 0.01_rt){
               return l1 + l3 * pow(x, 2); // Use Taylor expansion for small x
          }
          amrex::Real y =  exp(-2._rt * x);
          return (1._rt + y) / (1._rt - y) / x - 1._rt / pow(x, 2);
     };

    /**
     * \brief The derivative of the Langevin function.
     * 
     * \param x Input (unitless)
     * \return L'(x) = - 4 e^(-2x)/(1-e^(-2x)) + 1/x^2
     */
     static amrex::Real dLdx(amrex::Real x){
          if (x < 0.01_rt){
               return l1 + dldx2 * pow(x, 2); // Use Taylor expansion for small x
          }
          amrex::Real y =  exp(-2._rt * x);
          return -4._rt * y / pow(1._rt - y, 2) + 1._rt / pow(x, 2);
     };
};



#endif // WARPX_MACROSCOPIC_PROPERTIES_H_
