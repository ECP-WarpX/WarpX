#ifndef WARPX_JAMODEL_H_
#define WARPX_JAMODEL_H_

#include "MacroscopicProperties_fwd.H"

#include "Utils/WarpXConst.H"
#include "RealMatrix.H"

#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_Dim3.H>
#include <Dim3Util.H>

#include <memory>
#include <string>

/**
 * \brief Jiles-Atherton model of magnetic hysteresis.
 */
class JAModel
    : public MagneticMaterial
{
public:

    /**
     * \brief Specify the parameters of the Jiles-Atherton model for a magnetic material.
     *
     * \param alpha (unitless) Defines the effective field `H_e = H + alpha M`,
     *              where `H` is the magnetic field and `M` is the total magnetization.
     *              Quantifies interdomain coupling in the magnetic material.
     *
     * \param a (A/m) Describes the scale of the `H` field at which the magnetic material saturates
     *          via the ratio `|H_e|/a)`, where `H_e` is the effective field.
     *          Quantifies domain walls density in the magnetic material.
     *
     * \param Ms (A/m) The saturation magnetization of the magnetic material. The magnitude of
     *           the total magnetization `M` in the material should be bounded by `M_s`.
     *
     * \param k (A/m) Quantifies average energy required to break a pinning site in the magnetic material.
     *          Describes the scale by which the irreversible magnetization `M_irr` changes due to a change
     *          in the effective field `H_e`. In particular, `dM_irr/dH_e = (M_an - M_irr)/k`,
     *          where `M_an` is the anhysteretic magnetization.
     *
     * \param c (unitless) Quantifies magnetization reversibility: `M_rev = c (M_an - M_irr)`,
     *          where `M_rev`, `M_an`, and `M_irr` are the reversible, anhysteretic, and irreversible
     *          magnetizations, respectively. This equation corresponds to the bulge of the domain walls for
     *          small displacement. The total magnetization is `M = M_rev + M_irr = c M_an + (1-c) M_irr`.
     */
    JAModel (amrex::Real alpha,
             amrex::Real a,
             amrex::Real Ms,
             amrex::Real k,
             amrex::Real c);

    /**
     * \brief Update the magnetic field `H` and magnetization `M` for the Jiles-Atherton model.
     *
     * Both fields are needed as a way of keeping track of the previous state of the material.
     */
    void UpdateHandM (
              amrex::Real& H_x     ,       amrex::Real& H_y     ,        amrex::Real& H_z    ,
              amrex::Real& M_x     ,       amrex::Real& M_y     ,        amrex::Real& M_z    ,
        const amrex::Real& B_next_x, const amrex::Real& B_next_y, const amrex::Real& B_next_z) const override;


    //! (unitless) The `alpha` parameter of the Jiles-Atherton model.
    const amrex::Real& alpha () const { return m_param_alpha; }
    //! (A/m) The `a` parameter of the Jiles-Atherton model.
    const amrex::Real& a () const { return m_param_a; }
    //! (A/m) The `M_s` parameter of the Jiles-Atherton model.
    const amrex::Real& Ms () const { return m_param_Ms; }
    //! (A/m) The `k` parameter of the Jiles-Atherton model.
    const amrex::Real& k () const { return m_param_k; }
    //! (unitless) The `c` parameter of the Jiles-Atherton model.
    const amrex::Real& c () const { return m_param_c; }

protected:

    /**
     * \brief Compute the anhysteretic quantities `M_an` and `chi_an` for the JA model.
     */
    void ComputeAnhystereticContribution (
              amrex::Real& M_an_x,       amrex::Real& M_an_y,       amrex::Real& M_an_z,
              amrex::GpuArray<amrex::Real,6>& chi_an,
        const amrex::Real& H_e_x , const amrex::Real& H_e_y , const amrex::Real& H_e_z) const;

    /**
     * \brief Compute the change in the effective field `dH_e` from `dB`.
     *
     * Let `B` stand for `B/mu0` in the following:
     * The effective field is `H_e = H + alpha M = B - (1-alpha) M`.
     * Therefore `dB/dH_e = I + (1-alpha) chi`, where `chi = dM/dH_e`,
     * and thus `dH_e = (I + (1-alpha) chi)^{-1} dB`.
     */
    void Compute_dH_e (
              amrex::Real& dH_e_x    ,        amrex::Real& dH_e_y    ,         amrex::Real& dH_e_z,
        const amrex::GpuArray<amrex::Real,6>& dM_dH_e,
        const amrex::Real& dB_mu0inv_x, const amrex::Real& dB_mu0inv_y,  const amrex::Real& dB_mu0inv_z) const;

    const amrex::Real m_param_alpha;
    const amrex::Real m_param_a;
    const amrex::Real m_param_Ms;
    const amrex::Real m_param_k;
    const amrex::Real m_param_c;
    const amrex::Real m_param_one_minus_alpha;
    const amrex::Real m_param_one_minus_alpha_inv;
    const amrex::Real m_param_Ms_over_a;
    const amrex::Real m_param_c_Ms_over_a;
    const amrex::Real m_param_a_inv;
    const amrex::Real m_param_a2_inv;

};

namespace Langevin {

// public:

    using R = amrex::Real;

    /**
     * \brief The Langevin function `L(x) = 1/tanh(x) - 1/x`.
     *
     * Uses a polynomial approximation of `x/3 + (constant) x^3` for small x.
     * The constant is chosen to ensure continuity at the cutoff.
     */
    static amrex::Real L (amrex::Real x);

    /**
     * \brief `L(x)/x`, where `L(x) = 1/tanh(x) - 1/x` is the Langevin function.
     *
     * Uses a polynomial approximation of `1/3 + (constant) x^2` for small x.
     * The constant is chosen to ensure continuity at the cutoff.
     */
    static amrex::Real L_over_x (amrex::Real x);

    /**
     * \brief The derivative of the Langevin function: `L'(x) = -1/sinh^2(x) + 1/x^2`.
     *
     * Uses a polynomial approximation of `1/3 + (constant) x^2` for small x.
     * The constant is chosen to ensure continuity at the cutoff.
     */
    static amrex::Real dLdx (amrex::Real x);

    /**
     * \brief `(L'(x) - L(x)/x) / x^2`, where `L(x) = 1/tanh(x) - 1/x` is the Langevin function.
     *
     * Uses a polynomial approximation of `-2/45 + (constant) x^2` for small x.
     * The constant is chosen to ensure continuity at the cutoff.
     */
    static amrex::Real dLdx_minus_L_over_x_over_x2 (amrex::Real x);

// private:

    static constexpr double cutoff = 0.01;
    static constexpr double L_cutoff = 1.0/std::tanh(cutoff) - 1.0/cutoff;
    static constexpr double L_over_x_cutoff = L_cutoff/cutoff;
    static constexpr double dLdx_cutoff = -1.0/(std::sinh(cutoff)*std::sinh(cutoff)) + 1.0/(cutoff*cutoff);
    static constexpr double dLdx_minus_Loverx_over_x2_cutoff = (dLdx_cutoff - L_over_x_cutoff) / (cutoff*cutoff);

    // L(x) = 1/3 x - 1/45 x^3 + O(x^5)
    static constexpr R oneThird = R(1.0/3.0);

    //! Slightly different from -1/45 so that `L(x)` and `L_over_x(x)` are continuous at the cutoff.
    static constexpr R LangConst1 = R((L_over_x_cutoff - 1.0/3.0)/(cutoff*cutoff));

    //! Slightly different from -1/15 so that `dLdx(x)` is continuous at the cutoff.
    static constexpr R LangConst2 = R((dLdx_cutoff - 1.0/3.0)/(cutoff*cutoff));

    //! Slightly different from 8/945 so that `dLdx_minus_L_over_x_over_x2(x)` is continuous at the cutoff.
    static constexpr R LangConst3 = R((dLdx_minus_Loverx_over_x2_cutoff + 2.0/45.0)/(cutoff*cutoff));

};

#endif // WARPX_JAMODEL_H_
