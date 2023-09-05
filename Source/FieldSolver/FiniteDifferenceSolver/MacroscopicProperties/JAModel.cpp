#include <JAModel.H>

JAModelParameters::JAModelParameters (amrex::Real a_alpha, amrex::Real a_a, amrex::Real a_Ms, amrex::Real a_k, amrex::Real a_c)
{
    m_alpha = a_alpha;
    m_a = a_a;
    m_Ms = a_Ms;
    m_k = a_k;
    m_c = a_c;
    m_one_minus_alpha = 1._rt - a_alpha;
    m_Ms_over_a = a_Ms / a_a;
    m_c_times_Ms_over_a = a_c * a_Ms / a_a;
}

void
JAModel::UpdateM (amrex::RealVect& M, const amrex::RealVect& B, 
                const amrex::RealVect& dB, JAModelParameters* ja_model_parameters)
{
    /* 1 - alpha */
    amrex::Real oneMinusAlpha = ja_model_parameters->one_minus_alpha();

    /* The anhysteretic magnetization */
    amrex::RealVect M_an;
    /* Jacobian matrix chi_an = dM_an/dH_e. The anhysteretic part to chi */
    amrex::Array<amrex::Real, 6> chi_an;
    /* The effective field H_e = H + alpha M = mu0^{-1} B - (1-alpha) M */
    amrex::RealVect H_e = (1._rt / PhysConst::mu0) * B - oneMinusAlpha * M;
    /* Compute M_an and chi_an */ 
    JAModel::ComputeAnhystereticContribution(M_an, chi_an, H_e, ja_model_parameters);
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
        dH_e = Compute_dHe(chi_an, dBovermu0, oneMinusAlpha);
    }
    /* Update the magnetization */
    M = M + (1._rt / oneMinusAlpha) * (dBovermu0 - dH_e);
};

void
JAModel::ComputeAnhystereticContribution (amrex::RealVect& M_an,
                                     amrex::Array<amrex::Real, 6>& chi_an,
                                     const amrex::RealVect& H_e,
                                     JAModelParameters* ja_model_parameters)
{
    amrex::Real a = ja_model_parameters->a();
    /* c Ms / a */
    amrex::Real c_Ms_over_a = ja_model_parameters->c_times_Ms_over_a();
    amrex::Real H_e_mag = H_e.vectorLength(); // |H_e|
    if (H_e_mag <= amrex::Real(0.0)) { return; }
    amrex::Real x = H_e_mag / a; // x = |H_e| / a
    amrex::Real L_over_x = Langevin::L_over_x(x); // L(x)/x
    amrex::Real dLdx = Langevin::dLdx(x); // L'(x)
    M_an = ((ja_model_parameters->Ms_over_a()) * L_over_x) * H_e;
    amrex::Real delta = dLdx - L_over_x;  // delta = L'(x) - L(x)/x
    amrex::RealVect H_e_hat = H_e / H_e_mag; // Unit vector in direction of H_e
    // chi_an = dM_an / dH_e is a symmetric matrix, so only 6 components need to be specified.
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

amrex::RealVect
JAModel::Compute_dHe (amrex::Array<amrex::Real, 6> chi, amrex::RealVect dBovermu0, amrex::Real oneMinusAlpha)
{
    /* Jacobian matrix mu0^{-1} dB/dHe = I + (1 - alpha) chi */
    std::array<amrex::Real, 3UL> d;


    amrex::Array<amrex::Real, 6> chi_an;

    RealMatrix dBovermu0_dHe;

    amrex::Real offDiag;

    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        int jp1choose2 = j * (j + 1) / 2;
        dBovermu0_dHe[j][j] = 1 + oneMinusAlpha * chi[j + jp1choose2];
        for (int i = 0; i < j; ++i) {
            offDiag = oneMinusAlpha * chi[i + jp1choose2];
            dBovermu0_dHe[i][j] = offDiag;
            dBovermu0_dHe[j][i] = offDiag;
        }
    }
    /* Jacobian matrix mu0 dHe/dB */
    RealMatrix dHe_dBovermu0 = dBovermu0_dHe.Inverse();
    /* Could be more efficient by writing a function that solves x = A^{-1} b directly */
    return dHe_dBovermu0.Dot(dBovermu0);
}

amrex::Real
Langevin::L (amrex::Real x)
{
    if (x < 0.01_rt) { // Use Taylor expansion for small x
        return l1 * x + l3 * std::pow(x, 3); 
    }
    return 1._rt / std::tanh(x) - 1._rt / x;
}

amrex::Real
Langevin::L_over_x (amrex::Real x)
{
    if (x < 0.01_rt) { // Use Taylor expansion for small x
        return l1 + l3 * x * x; 
    }
    amrex::Real xinv = 1._rt / x;
    return xinv * (1._rt / std::tanh(x) - xinv );
};

amrex::Real
Langevin::dLdx (amrex::Real x)
{
    if (x < 0.01_rt) { // Use Taylor expansion for small x
        return l1 + dldx2 * x * x; 
    }
    amrex::Real y = std::sinh(x);
    return -1._rt / (y * y) + 1._rt / (x * x);
};

