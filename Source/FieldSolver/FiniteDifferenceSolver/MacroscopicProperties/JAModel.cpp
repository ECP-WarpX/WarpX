#include <JAModel.H>

JAModel::JAModel (amrex::Real alpha,
                  amrex::Real a,
                  amrex::Real Ms,
                  amrex::Real k,
                  amrex::Real c)
{
    m_param_alpha = alpha;
    m_param_a = a;
    m_param_Ms = Ms;
    m_param_k = k;
    m_param_c = c;
    m_param_one_minus_alpha = amrex::Real(1.0) - alpha;
    m_param_one_minus_alpha_inv = amrex::Real(1.0) / m_param_one_minus_alpha;
    m_param_Ms_over_a = Ms / a;
    m_param_c_Ms_over_a = c * Ms / a;
    m_param_a_inv = amrex::Real(1.0) / a;
    m_param_a2_inv = m_param_a_inv*m_param_a_inv;
}

void JAModel::UpdateHandM (
          amrex::Real& H_x     ,       amrex::Real& H_y     ,       amrex::Real& H_z     ,
          amrex::Real& M_x     ,       amrex::Real& M_y     ,       amrex::Real& M_z     ,
    const amrex::Real& B_next_x, const amrex::Real& B_next_y, const amrex::Real& B_next_z) const
{
    const amrex::Real B_prev_mu0inv_x = H_x + M_x;
    const amrex::Real B_prev_mu0inv_y = H_y + M_y;
    const amrex::Real B_prev_mu0inv_z = H_z + M_z;

    const amrex::Real B_next_mu0inv_x = PhysConst::mu0inv * B_next_x;
    const amrex::Real B_next_mu0inv_y = PhysConst::mu0inv * B_next_y;
    const amrex::Real B_next_mu0inv_z = PhysConst::mu0inv * B_next_z;

    const amrex::Real dB_mu0inv_x = B_next_mu0inv_x - B_prev_mu0inv_x;
    const amrex::Real dB_mu0inv_y = B_next_mu0inv_y - B_prev_mu0inv_y;
    const amrex::Real dB_mu0inv_z = B_next_mu0inv_z - B_prev_mu0inv_z;

    /* The effective field is H_e = H + alpha M = mu0^{-1} B - (1-alpha) M */
    amrex::Real H_e_x = B_prev_mu0inv_x - m_param_one_minus_alpha * M_x;
    amrex::Real H_e_y = B_prev_mu0inv_y - m_param_one_minus_alpha * M_y;
    amrex::Real H_e_z = B_prev_mu0inv_z - m_param_one_minus_alpha * M_z;

    /* The anhysteretic magnetization */
    amrex::Real M_an_x, M_an_y, M_an_z;

    /* Jacobian matrix chi_an = dM_an/dH_e. The anhysteretic contribution to chi */
    amrex::GpuArray<amrex::Real,6> chi_an;

    /* Compute M_an and chi_an */
    JAModel::ComputeAnhystereticContribution(M_an_x,M_an_y,M_an_z,chi_an,H_e_x,H_e_y,H_e_z);

    //! M_diff = M_an - M
    const amrex::Real M_diff_x = M_an_x - M_x;
    const amrex::Real M_diff_y = M_an_y - M_y;
    const amrex::Real M_diff_z = M_an_z - M_z;

    //! |M_an - M|
    const amrex::Real M_diff_mag = amrex::Real(std::sqrt(M_diff_x*M_diff_x + M_diff_y*M_diff_y + M_diff_z*M_diff_z));

    const amrex::Real denominator = m_param_k * M_diff_mag; // k |M_an - M|

    //! The change in the effective field
    amrex::Real dH_e_x, dH_e_y, dH_e_z;

    if (denominator != amrex::Real(0.0)) {
        //! Jacobian matrix chi = dM/dH_e = chi_an + chi_irr. Is symmetric 3-by-3 so only need 6 components
        amrex::GpuArray<amrex::Real,6> chi;
        chi[0] = chi_an[0] + (M_diff_x*M_diff_x / denominator); // xx
        chi[1] = chi_an[0] + (M_diff_x*M_diff_y / denominator); // xy
        chi[2] = chi_an[0] + (M_diff_y*M_diff_y / denominator); // yy
        chi[3] = chi_an[3] + (M_diff_x*M_diff_z / denominator); // xz
        chi[4] = chi_an[4] + (M_diff_y*M_diff_z / denominator); // yz
        chi[5] = chi_an[5] + (M_diff_z*M_diff_z / denominator); // zz
        /* Compute dH_e using the total contribution */
        Compute_dH_e(dH_e_x,dH_e_y,dH_e_z,chi,dB_mu0inv_x,dB_mu0inv_y,dB_mu0inv_z);
        /* Check the sign of (M_an - M) dot dH_e */
        if (M_diff_x*dH_e_x + M_diff_y*dH_e_y + M_diff_z*dH_e_z < Real(0.0)) {
            /* Use only anhysteretic part chi_an if (M_an - M) dot dH_e < 0 */
            Compute_dH_e(dH_e_x,dH_e_y,dH_e_z,chi_an,dB_mu0inv_x,dB_mu0inv_y,dB_mu0inv_z);
        }
    }
    else {
        Compute_dH_e(dH_e_x,dH_e_y,dH_e_z,chi_an,dB_mu0inv_x,dB_mu0inv_y,dB_mu0inv_z);
    }

    //! Update the magnetization. dM = (dB/mu0 - dH_e)/(1-alpha)
    M_x += m_param_one_minus_alpha_inv * (dB_mu0inv_x - dH_e_x);
    M_y += m_param_one_minus_alpha_inv * (dB_mu0inv_y - dH_e_y);
    M_z += m_param_one_minus_alpha_inv * (dB_mu0inv_z - dH_e_z);

    //! Update the magnetic field
    H_x = B_next_mu0inv_x - M_x;
    H_y = B_next_mu0inv_y - M_y;
    H_z = B_next_mu0inv_z - M_z;

}

void JAModel::ComputeAnhystereticContribution (
          amrex::Real& M_an_x,       amrex::Real& M_an_y,       amrex::Real& M_an_z,
          amrex::GpuArray<amrex::Real,6>& chi_an,
    const amrex::Real& H_e_x , const amrex::Real& H_e_y , const amrex::Real& H_e_z ) const
{
    //! |H_e|
    amrex::Real H_e_mag = amrex::Real(std::sqrt(H_e_x*H_e_x + H_e_y*H_e_y + H_e_z*H_e_z));
    //! |H_e| / a
    amrex::Real x = H_e_mag / m_param_a;
    //! L(x) / x, where x = |H_e| / a
    amrex::Real L_over_x = Langevin::L_over_x(x);
    //! (L'(x) - L(x)/x) / (a^2 x^2) , where x = |H_e| / a
    amrex::Real dLdx_minus_L_over_x_over_x2a2 = Langevin::dLdx_minus_L_over_x_over_x2(x) / (m_param_a*m_param_a);

    //! M_an = M_s L(|H_e|/a) H_e_hat. It is written with L_over_x to avoid errors if H_e_mag is too small
    M_an_x = m_param_Ms_over_a * L_over_x * H_e_x;
    M_an_y = m_param_Ms_over_a * L_over_x * H_e_y;
    M_an_z = m_param_Ms_over_a * L_over_x * H_e_z;

    // chi_an = dM_an / dH_e is a symmetric matrix, so only 6 components need to be specified.
    // The order of elements is xx, xy, yy, xz, yz, zz)
    chi_an[0] = m_param_c_Ms_over_a * (L_over_x + dLdx_minus_L_over_x_over_x2a2 * H_e_x * H_e_x);
    chi_an[1] = m_param_c_Ms_over_a * (           dLdx_minus_L_over_x_over_x2a2 * H_e_x * H_e_y);
    chi_an[2] = m_param_c_Ms_over_a * (L_over_x + dLdx_minus_L_over_x_over_x2a2 * H_e_y * H_e_y);
    chi_an[3] = m_param_c_Ms_over_a * (           dLdx_minus_L_over_x_over_x2a2 * H_e_x * H_e_z);
    chi_an[4] = m_param_c_Ms_over_a * (           dLdx_minus_L_over_x_over_x2a2 * H_e_y * H_e_z);
    chi_an[5] = m_param_c_Ms_over_a * (L_over_x + dLdx_minus_L_over_x_over_x2a2 * H_e_z * H_e_z);
}

void JAModel::Compute_dH_e (
          amrex::Real& dH_e_x    ,        amrex::Real& dH_e_y    ,         amrex::Real& dH_e_z     ,
    const amrex::GpuArray<amrex::Real,6>& chi_an,
    const amrex::Real& dB_mu0inv_x, const amrex::Real& dB_mu0inv_y,  const amrex::Real& dB_mu0inv_z) const
{
    /* Jacobian matrix mu0^{-1} dB/dHe = I + (1-alpha) chi */
    const amrex::Real dBdHe00 = 1 + m_param_one_minus_alpha * chi[0];
    const amrex::Real dBdHe01 =     m_param_one_minus_alpha * chi[1];
    const amrex::Real dBdHe11 = 1 + m_param_one_minus_alpha * chi[2];
    const amrex::Real dBdHe02 =     m_param_one_minus_alpha * chi[3];
    const amrex::Real dBdHe12 =     m_param_one_minus_alpha * chi[4];
    const amrex::Real dBdHe22 = 1 + m_param_one_minus_alpha * chi[5];

    const amrex::Real det = dBdHe00*dBdHe11*dBdHe22 + amrex::Real(2.0)*dBdHe01*dBdHe02*dBdHe12
                          - dBdHe12*dBdHe12*dBdHe00 - dBdHe02*dBdHe02*dBdHe11 - dBdHe01*dBdHe01*dBdHe22;
    const amrex::Real detinv = amrex::Real(1.0) / det;

    const amrex::Real dHedB00 = (dBdHe11*dBdHe22 - dBdHe12*dBdHe12)*detinv;
    const amrex::Real dHedB01 = (dBdHe02*dBdHe12 - dBdHe01*dBdHe22)*detinv;
    const amrex::Real dHedB11 = (dBdHe00*dBdHe22 - dBdHe02*dBdHe02)*detinv;
    const amrex::Real dHedB02 = (dBdHe01*dBdHe12 - dBdHe02*dBdHe11)*detinv;
    const amrex::Real dHedB12 = (dBdHe01*dBdHe02 - dBdHe00*dBdHe12)*detinv;
    const amrex::Real dHedB22 = (dBdHe00*dBdHe11 - dBdHe01*dBdHe01)*detinv;

    dH_e_x = dHedB00*dB_mu0inv_x + dHedB01*dB_mu0inv_y + dHedB02*dB_mu0inv_z;
    dH_e_y = dHedB01*dB_mu0inv_x + dHedB11*dB_mu0inv_y + dHedB12*dB_mu0inv_z;
    dH_e_z = dHedB02*dB_mu0inv_x + dHedB12*dB_mu0inv_y + dHedB22*dB_mu0inv_z;
}

amrex::Real
Langevin::L (amrex::Real x)
{
    if (x < R(cutoff)) { // Use Taylor expansion for small x
        return (oneThird + LangConst1*x*x)*x;
    }
    return R(1.0 / std::tanh(x)) - R(1.0) / x;
}

amrex::Real
Langevin::L_over_x (amrex::Real x)
{
    if (x < R(cutoff)) { // Use Taylor expansion for small x
        return oneThird + LangConst1*x*x;
    }
    R xinv = R(1.0) / x;
    return xinv*(R(1.0 / std::tanh(x)) - xinv);
};

amrex::Real
Langevin::dLdx (amrex::Real x)
{
    if (x < R(cutoff)) { // Use Taylor expansion for small x
        return oneThird + LangConst2 * x * x;
    }
    const auto sinh = std::sinh(x);
    return R(-1.0 / (sinh * sinh)) + R(1.0) / (x * x);
};

amrex::Real
Langevin::dLdx_minus_L_over_x_over_x2 (amrex::Real x)
{
    if (x < R(cutoff)) {
        return R(-2.0/45.0) + LangConst3*x*x;
    }
    const R xinv = R(1.0) / x;
    const R xinv2 = xinv*xinv;
    const auto csch = 1.0 / std::sinh(x);
    const auto coth = 1.0 / std::tanh(x);
    return xinv2*(R(2.0)*xinv2 - xinv*R(coth) - R(csch*csch));
}
