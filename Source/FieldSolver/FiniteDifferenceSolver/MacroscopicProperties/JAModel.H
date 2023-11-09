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
 * \brief This struct contains the five parameters necessary to describe magnetic
 * hysteresis of a material using the Jiles-Atherton (JA) model.
 */
struct JAModelParameters {

public:
    /**
     * \brief Specify the parameters of the Jiles-Atherton (JA) model for a magnetic material.
     *
     * \param alpha (unitless) Defines the effective field `H_e = H + alpha M`,
     *              where `H` is the magnetic field and `M` is the total magnetization.
     *              Quantifies interdomain coupling in the magnetic material.
     *
     * \param a (A/m) Describes the scale of the `H` field at which the magnetic
     *          material saturates. Used in the Langevin function as `L(H_e / a)`.
     *          Quantifies domain walls density in the magnetic material.
     *
     * \param Ms (A/m) Saturation magnetization of the magnetic material.
     *
     * \param k (A/m) Quantifies average energy required to break a pinning site in the magnetic material.
     *          Describes the scale by which the irreversible magnetization `M_irr` changes due to a change
     *          in the effective field `H_e`. In particular, `dM_irr / dH_e = (M_an - M_irr) / k`,
     *          where `M_an` is the anhysteretic magnetization.
     *
     * \param c (unitless) Quantifies magnetization reversibility: `M_rev = c (M_an - M_irr)`,
     *          where `M_rev`, `M_an`, and `M_irr` are the reversible, anhysteretic, and irreversible
     *          magnetizations, respectively. This equation corresponds to the bulge of the domain walls
     *          for small displacement. The total magnetization is
     *          `M = M_rev + M_irr = c M_an + (1-c) M_irr`.
     */
    JAModelParameters (amrex::Real a_alpha, amrex::Real a_a, amrex::Real a_Ms, amrex::Real a_k, amrex::Real a_c);

    /**
     * \brief The alpha parameter in the Jiles-Atherton Model
     *
     * (unitless) Defines the effective field `H_e = H + alpha M`,
     * where `H` is the magnetic field and `M` is the total magnetization.
     * Quantifies interdomain coupling in the magnetic material.
     */
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

class JAModel
    : public MagneticMaterial
{
public:

    using GpuArray3 = amrex::GpuArray<amrex::Real,3>;
    using GpuArray6 = amrex::GpuArray<amrex::Real,6>;
    using Vec3 = amrex::Dim3;
    using RT = amrex::Real;

    /**
     * \brief Update the magnetization and magnetic field for the Jiles-Atherton model
     */
    void
    UpdateHandM (      RT& H_x    ,       RT& H_y    ,       RT& H_z    ,
                       RT& M_x    ,       RT& M_y    ,       RT& M_z    ,
                 const RT& Bnext_x, const RT& Bnext_y, const RT& Bnext_z) const override;

    // void ComputeH(GpuArray3 &H, const GpuArray3 &B, const GpuArray3 &dB);
    // void ComputeH(RT &Hx, RT &Hy, RT &Hz,
    //               const RT &Bx, const RT &By, const RT &Bz,
    //               const RT &dBx, const GpuArray3 &dB);

    // void ComputeH (      RT& H_x    ,       RT& H_y    ,       RT& H_z    ,
    //                      RT& M_x    ,       RT& M_y    ,       RT& M_z    ,
    //                const RT& Bprev_muin0_x, const RT& Bprev_muin0_y, const RT& Bprev_muin0_z,
    //                const RT& Bnext_muin0_x, const RT& Bnext_muin0_y, const RT& Bnext_muin0_z) const override;

    // void UpdateM (      RT& M_x          ,       RT& M_y          ,       RT& M_z          ,
    //               const RT& Bprev_muin0_x, const RT& Bprev_muin0_y, const RT& Bprev_muin0_z,
    //               const RT& Bnext_muin0_x, const RT& Bnext_muin0_y, const RT& Bnext_muin0_z) const override;

    // void UpdateMx (      RT& Mx,  const RT& My,  const RT& Mz,
    //                const RT& Bx,  const RT& By,  const RT& Bz,
    //                const RT& dBx, const RT& dBy, const RT& dBz);

protected:

    /**
     * \brief Compute the anhysteretic quantities for the JA model.
     *
     * \param[out] M_an The anhysteretic magnetization
     * \param[out] chi_an The Jacobian matrix dM_an/dH_e. Symmetric, with elements ordered as xx, xy, yy, xz, yz, zz.
     * \param[in] H_e The effective magnetic field
     */
    void ComputeAnhystereticContribution (      RT& M_an_x,       RT& M_an_y,       RT& M_an_z,
                                          amrex::GpuArray<RT,6>& chi_an,
                                          const RT& H_e_x , const RT& H_e_y , const RT& H_e_z) const;

    /**
     * \brief Computes change in effective field: dH_e = (dB/dH_e)^{-1} dB.
     *
     * Since B/mu0 = H_e + (1-alpha) M, we have d(B/mu0)/dH_e = I + (1-alpha) chi,
     * where chi = dM/dH_e.
     */
    Compute_dHe (      GpuArray3 &dH_e,
                 const GpuArray6 &chi,
                 const GpuArray3 &dB_over_mu0);

    JAModelParameters m_ja_model_parameters;


};

struct Langevin {

public:

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

private:

    static constexpr double cutoff = 0.01;
    static constexpr double L_cutoff = 1.0/std::tanh(cutoff) - 1.0/cutoff;
    static constexpr double L_over_x_cutoff = L_cutoff / cutoff;
    static constexpr double dLdx_cutoff = -1.0/(std::sinh(cutoff)*std::sinh(cutoff)) + 1.0/(cutoff*cutoff);
    static constexpr double dLdx_minus_Loverx_over_x2_cutoff = (dLdx_cutoff - L_over_x_cutoff) / (cutoff*cutoff);

    // L(x) = 1/3 x - 1/45 x^3 + O(x^5)
    static constexpr R oneThird = R(1.0/3.0);

    //! Slightly different from -1/45 so that `L(x)` and `L_over_x(x)` are continuous at the cutoff.
    static constexpr R LangConst1 = R((L_over_x_cutoff - 1.0/3.0) / (cutoff*cutoff));

    //! Slightly different from -1/15 so that `dLdx(x)` is continuous at the cutoff.
    static constexpr R LangConst2 = R((dLdx_cutoff - 1.0/3.0) / (cutoff*cutoff));

    //! Slightly different from 8/945 so that `dLdx_minus_L_over_x_over_x2(x)` is continuous at the cutoff.
    static constexpr R LangConst3 = R((dLdx_minus_Loverx_over_x2_cutoff + 2.0/45.0) / (cutoff*cutoff));


    // static constexpr amrex::Real LangConst1 = -0.02222201058412696_rt;
    // static constexpr double LangConst1 = -0.02222201058412696;
    // static constexpr double LangConst1 = -2.2222201058412696E-2;
    // static constexpr double LangConst1 = (L_over_x_cutoff - 1.0/3.0) / (cutoff*cutoff);
    // static constexpr amrex::Real LangConst2 = -0.06666559177859188_rt;
    // static constexpr double LangConst2 = -0.06666560848042309;
    // static constexpr double LangConst2 = -6.666560848042309E-2;
    // static constexpr double LangConst2 = (dLdx_cutoff - 1.0/3.0) / (cutoff*cutoff);
    // static constexpr amrex::Real LangConst3 = 0.008665953929384851;
    // static constexpr double LangConst3 = 8.0/945.0 - (2.0/1575.0)*cutoff*cutoff + (16.0/93555.0)*cutoff*cutoff*cutoff*cutoff;
    // static constexpr double LangConst3 = 0.00846548148319168;
    // static constexpr double LangConst3 = 8.46548148319168E-3;
    // static constexpr double LangConst3 = (dLdx_minus_Loverx_over_x2_cutoff + 2.0/45.0) / (cutoff*cutoff);

};

