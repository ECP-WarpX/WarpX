#ifndef WARPX_JAMODEL_H_
#define WARPX_JAMODEL_H_

#include "MacroscopicProperties_fwd.H"

#include "Utils/WarpXConst.H"

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>

/**
 * \brief Jiles-Atherton model of magnetic hysteresis.
 */
// : public MagneticMaterial
class JAModel {
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

    amrex::Real m_param_alpha;
    amrex::Real m_param_a;
    amrex::Real m_param_Ms;
    amrex::Real m_param_k;
    amrex::Real m_param_c;

    amrex::Real m_param_one_minus_alpha;
    amrex::Real m_param_one_minus_alpha_inv;
    amrex::Real m_param_Ms_over_a;
    amrex::Real m_param_c_Ms_over_a;
    amrex::Real m_param_a_inv;
    amrex::Real m_param_a2_inv;

};

// M330-50 A steel Jiles-Atherton model parameters. Source:
// Chwastek K. The applications of fixed‐point theorem in optimisation problems.
// Arch Electr Eng. 2012;61(2):189–198
// static const JAModel mag_mat_M330_50_A_steel{amrex::Real(7.45E-5) /* alpha (unitless) */,
//                                              amrex::Real(26.1)    /* a (A/m) */,
//                                              amrex::Real(1.26E6)  /* Ms (A/m) */,
//                                              amrex::Real(52.3)    /* k (A/m) */,
//                                              amrex::Real(0.13)    /* c (unitless) */};

#endif // WARPX_JAMODEL_H_
