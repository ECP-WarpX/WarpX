#include "MacroscopicProperties_fwd.H"

#include "Utils/WarpXConst.H"
#include "RealMatrix.H"

#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>

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
    JAModelParameters (amrex::Real a_alpha, amrex::Real a_a, amrex::Real a_Ms, amrex::Real a_k, amrex::Real a_c);

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

    /**(unitless) Quantifies magnetization reversiblity: `M_rev = c (M_an - M_irr)`,
     * where `M_rev`, `M_an`, and `M_irr` are the reversible, anhysteretic, and irreversible
     * magnetizations, respectively. This equation corresponds to the bulge of the domain walls 
     * for small displacement. The total magnetization is 
     * `M = M_rev + M_irr = c M_an + (1-c) M_irr`.
     */
    amrex::Real c () { return m_c; }

    amrex::Real one_minus_alpha () { return m_one_minus_alpha; }

    amrex::Real Ms_over_a () { return m_Ms_over_a; }

    amrex::Real c_times_Ms_over_a () { return m_c_times_Ms_over_a; }

private:

    amrex::Real m_alpha;
    amrex::Real m_a;
    amrex::Real m_Ms;
    amrex::Real m_k;
    amrex::Real m_c;
    amrex::Real m_one_minus_alpha;
    amrex::Real m_Ms_over_a;
    amrex::Real m_c_times_Ms_over_a;

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
            const amrex::RealVect &dB, JAModelParameters *ja_model_parameters);

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
                                     JAModelParameters* ja_model_parameters);

    /**
     * \brief
     * Computes change in effective field as dH_e = (dB/dH_e)^{-1} dB.
     * Since B/mu0 = H_e + (1-alpha) M, we have d(B/mu0)/dH_e = I + (1-alpha) chi, 
     * where chi = dM/dH_e.
     * 
     * \param chi Components of symmetric matrix chi = dM/dH_e
     * \param dB_over_mu0 mu0^{-1} dB
     * \param oneMinusAlpha 1 - alpha
     * \return amrex::RealVect dH_e = (dB/dH_e)^{-1} dB
     */
    amrex::RealVect
    Compute_dHe (amrex::Array<amrex::Real, 6> chi, amrex::RealVect dB_over_mu0, amrex::Real oneMinusAlpha);

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
    L (amrex::Real x);

    /**
     * \brief The Langevin function L(x) over x.
     * 
     * \param x Input (unitless)
     * \return L(x)/x = coth(x)/x - 1/x^2
     */
     static amrex::Real
     L_over_x (amrex::Real x);

    /**
     * \brief The derivative of the Langevin function.
     * 
     * \param x Input (unitless)
     * \return L'(x) = - 1/sinh^2(x) + 1/x^2
     */
     static amrex::Real
     dLdx (amrex::Real x);
};

