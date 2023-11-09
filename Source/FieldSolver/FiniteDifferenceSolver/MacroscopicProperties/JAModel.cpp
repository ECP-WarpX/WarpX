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
JAModel::UpdateHandM (      RT& H_x    ,       RT& H_y    ,       RT& H_z    ,
                            RT& M_x    ,       RT& M_y    ,       RT& M_z    ,
                      const RT& Bnext_x, const RT& Bnext_y, const RT& Bnext_z) const
{
    const RT Bprev_mu0inv_x = H_x + M_x;
    const RT Bprev_mu0inv_y = H_y + M_y;
    const RT Bprev_mu0inv_z = H_z + M_z;

    const RT Bnext_mu0inv_x = PhysConst::mu0inv * Bnext_x;
    const RT Bnext_mu0inv_y = PhysConst::mu0inv * Bnext_y;
    const RT Bnext_mu0inv_z = PhysConst::mu0inv * Bnext_z;

    const RT dB_mu0inv_x = Bnext_mu0inv_x - Bprev_mu0inv_x;
    const RT dB_mu0inv_y = Bnext_mu0inv_y - Bprev_mu0inv_y;
    const RT dB_mu0inv_z = Bnext_mu0inv_z - Bprev_mu0inv_z;

    /* The effective field is H_e = H + alpha M = mu0^{-1} B - (1-alpha) M */
    RT H_e_x = Bprev_mu0inv_x - m_param_one_minus_alpha * M_x;
    RT H_e_y = Bprev_mu0inv_y - m_param_one_minus_alpha * M_y;
    RT H_e_z = Bprev_mu0inv_z - m_param_one_minus_alpha * M_z;

    /* The anhysteretic magnetization */
    RT M_an_x, M_an_y, M_an_z;

    /* Jacobian matrix chi_an = dM_an/dH_e. The anhysteretic contribution to chi */
    GpuArray6 chi_an;

    /* Compute M_an and chi_an */
    JAModel::ComputeAnhystereticContribution(M_an_x,M_an_y,M_an_z,chi_an,H_e_x,H_e_y,H_e_z);

    //! M_diff = M_an - M
    const RT M_diff_x = M_an_x - M_x;
    const RT M_diff_y = M_an_y - M_y;
    const RT M_diff_z = M_an_z - M_z;

    //! |M_an - M|
    const RT M_diff_mag = RT(std::sqrt(M_diff_x*M_diff_x + M_diff_y*M_diff_y + M_diff_z*M_diff_z));

    const RT denominator = m_param_k * M_diff_mag; // k |M_an - M|

    //! The change in the effective field
    RT dH_e_x, dH_e_y, dH_e_z;

    if (denominator != RT(0.0)) {
        //! Jacobian matrix chi = dM/dH_e = chi_an + chi_irr. Is symmetric 3-by-3 so only need 6 components
        GpuArray6 chi;
        chi[0] = chi_an[0] + (M_diff_x*M_diff_x / denominator); // xx
        chi[1] = chi_an[0] + (M_diff_x*M_diff_y / denominator); // xy
        chi[2] = chi_an[0] + (M_diff_y*M_diff_y / denominator); // yy
        chi[3] = chi_an[3] + (M_diff_x*M_diff_z / denominator); // xz
        chi[4] = chi_an[4] + (M_diff_y*M_diff_z / denominator); // yz
        chi[5] = chi_an[5] + (M_diff_z*M_diff_z / denominator); // zz
// #endif
        /* Compute dH_e using the total contribution */
        Compute_dHe(dH_e_x,dH_e_y,dH_e_z,chi,dB_mu0inv_x,dB_mu0inv_y,dB_mu0inv_z);
        /* Check the sign of (M_an - M) dot dH_e */
        if (M_diff_x*dH_e_x + M_diff_y*dH_e_y + M_diff_z*dH_e_z < RT(0.0)) {
            /* Use only anhysteretic part chi_an if (M_an - M) dot dH_e < 0 */
            Compute_dHe(dH_e_x,dH_e_y,dH_e_z,chi_an,dB_mu0inv_x,dB_mu0inv_y,dB_mu0inv_z);
        }
    }
    else {
        Compute_dHe(dH_e_x,dH_e_y,dH_e_z,chi_an,dB_mu0inv_x,dB_mu0inv_y,dB_mu0inv_z);
    }

    //! Update the magnetization. dM = (dB/mu0 - dH_e)/(1-alpha)
    M_x += m_param_one_minus_alpha_inv * (dB_mu0inv_x - dH_e_x);
    M_y += m_param_one_minus_alpha_inv * (dB_mu0inv_y - dH_e_y);
    M_z += m_param_one_minus_alpha_inv * (dB_mu0inv_z - dH_e_z);
}

void
JAModel::ComputeAnhystereticContribution (      RT& M_an_x,       RT& M_an_y,       RT& M_an_z,
                                          amrex::GpuArray<RT,6>& chi_an,
                                          const RT& H_e_x , const RT& H_e_y , const RT& H_e_z) const
{
    using RT = amrex::Real;
    //! |H_e|
    RT H_e_mag = RT(std::sqrt(H_e_x*H_e_x + H_e_y*H_e_y + H_e_z*H_e_z));
    //! |H_e| / a
    RT x = H_e_mag / m_param_a;
    //! L(x) / x, where x = |H_e| / a
    RT L_over_x = Langevin::L_over_x(x);
    //! (L'(x) - L(x)/x) / (a^2 x^2) , where x = |H_e| / a
    RT delta = Langevin::dLdx_minus_L_over_x_over_x2(x) / (m_param_a*m_param_a);

    //! M_an = M_s L(|H_e|/a) H_e_hat. It is written with L_over_x to avoid errors if H_e_mag is too small
    M_an_x = m_param_Ms_over_a * L_over_x * H_e_x;
    M_an_y = m_param_Ms_over_a * L_over_x * H_e_y;
    M_an_z = m_param_Ms_over_a * L_over_x * H_e_z;

    // chi_an = dM_an / dH_e is a symmetric matrix, so only 6 components need to be specified.
    // The order of elements is xx, xy, yy, xz, yz, zz)
    chi_an[0] = m_param_c_Ms_over_a * (L_over_x + delta * H_e_x * H_e_x);
    chi_an[1] = m_param_c_Ms_over_a * (           delta * H_e_x * H_e_y);
    chi_an[2] = m_param_c_Ms_over_a * (L_over_x + delta * H_e_y * H_e_y);
    chi_an[3] = m_param_c_Ms_over_a * (           delta * H_e_x * H_e_z);
    chi_an[4] = m_param_c_Ms_over_a * (           delta * H_e_y * H_e_z);
    chi_an[5] = m_param_c_Ms_over_a * (L_over_x + delta * H_e_z * H_e_z);
}

void
JAModel::Compute_dHe (      RT& dH_e_x    ,        RT& dH_e_y    ,         RT& dH_e_z     ,
                      const amrex::GpuArray<RT,6>& chi_an                                 ,
                      const RT& dB_mu0inv_x, const RT& dB_mu0inv_y,  const RT& dB_mu0inv_z)
{
    /* Jacobian matrix mu0^{-1} dB/dHe = I + (1 - alpha) chi */
    RealMatrix dBovermu0_dHe;

    amrex::Real dBovermu0_dHe_00 = 1 + m_param_oneMinusAlpha * chi[0];
    amrex::Real dBovermu0_dHe_01 =     oneMinusAlpha * chi[1];
    amrex::Real dBovermu0_dHe_11 = 1 + oneMinusAlpha * chi[2];
    amrex::Real dBovermu0_dHe_02 =     oneMinusAlpha * chi[3];
    amrex::Real dBovermu0_dHe_12 =     oneMinusAlpha * chi[4];
    amrex::Real dBovermu0_dHe_22 = 1 + oneMinusAlpha * chi[5];

    amrex::Real factor = dBovermu0_dHe_12 / dBovermu0_dHe_11;

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

// void
// JAModel::ComputeH (      RT& H_x    ,       RT& H_y    ,       RT& H_z    ,
//                          RT& M_x    ,       RT& M_y    ,       RT& M_z    ,
//                    const RT& Bprev_x, const RT& Bprev_y, const RT& Bprev_z,
//                    const RT& Bnext_x, const RT& Bnext_y, const RT& Bnext_z) const
// {
//     Vec3 H = {H_x,H_y,H_z};
//     UpdateM(M_x,M_y,M_z,Bprev_x,Bprev_y,Bprev_z,Bnext_x,Bnext_y,Bnext_z);
//     H_x = PhysConst::mu0inv * Bprev_x
// }

// void
// JAModel::UpdateM (      RT& M_x ,       RT &M_y,        RT &M_z,
//                   const RT& B_x , const RT &B_y , const RT &B_z,
//                   const RT& dB_x, const RT &dB_y, const RT &dB_z)
// {
//     /* 1 - alpha */
//     amrex::Real oneMinusAlpha = m_ja_model_parameters.one_minus_alpha();

//     /* The anhysteretic magnetization */
//     GpuArray3 M_an;
//     /* Jacobian matrix chi_an = dM_an/dH_e. The anhysteretic contribution to chi */
//     GpuArray6 chi_an;
//     /* The effective field H_e = H + alpha M = mu0^{-1} B - (1-alpha) M */
//     GpuArray3 H_e;
//     H_e[0] = PhysConst::mu0inv * B_x - oneMinusAlpha * M_x;
//     H_e[1] = PhysConst::mu0inv * B_y - oneMinusAlpha * M_y;
//     H_e[2] = PhysConst::mu0inv * B_z - oneMinusAlpha * M_z;

//     /* Compute M_an and chi_an */
//     JAModel::ComputeAnhystereticContribution(M_an, chi_an, H_e);

//     /* The change in the effective field */
//     GpuArray3 dH_e;
//     GpuArray3 dBovermu0 = PhysConst::mu0inv * dB; // mu0^{-1} dB
//     GpuArray3 M_diff = M_an - M; // M_an - M


//     amrex::Real M_diff_mag = M_diff.vectorLength(); // |M_an - M|
//     amrex::Real denominator = m_ja_model_parameters.k() * M_diff_mag; // k |M_an - M|

//     if (denominator != 0) {
//         // Jacobian matrix chi = dM/dH_e = chi_an + chi_irr
//         amrex::Array<amrex::Real, 6> chi(chi_an);
//         // Add the irreversible contribution chi_irr.
//         chi[0] += M_diff[0] * M_diff[0] / denominator; // xx
// // #if AMREX_SPACEDIM > 1
//         chi[1] += M_diff[1] * M_diff[0] / denominator; // xy
//         chi[2] += M_diff[1] * M_diff[1] / denominator; // yy
// // #endif
// // #if AMREX_SPACEDIM > 2
//         chi[3] += M_diff[2] * M_diff[0] / denominator; // xz
//         chi[4] += M_diff[2] * M_diff[1] / denominator; // yz
//         chi[5] += M_diff[2] * M_diff[2] / denominator; // zz
// // #endif
//         /* Compute dH_e using the total contribution */
//         Compute_dHe(dHe, chi, dBovermu0, oneMinusAlpha);
//         /* Check the sign of (M_an - M) dot dH_e */
//         if (M_diff.dotProduct(dH_e) < 0) {
//             /* Use only anhysteretic part of chi if (M_an - M) dot dH_e < 0 */
//             Compute_dHe(dHe, chi_an, dBovermu0, oneMinusAlpha);
//         }
//     }
//     else {
//         Compute_dHe(dHe, chi_an, dBovermu0, oneMinusAlpha);
//     }

//     /* Update the magnetization */
//     M = M + (1._rt / oneMinusAlpha) * (dBovermu0 - dH_e);


//     M = amrex::Dim3_plus(M,amrex::Dim3_mult(1._rt / oneMinusAlpha, amrex::Dim3_minus(dBovermu0, dH_e)));
// }

// void
// JAModel::UpdateMx (      RT& M_x,  const RT& M_y,  const RT& M_z,
//                    const RT& Bprev_x,  const RT& Bprev_y,  const RT& Bprev_z,
//                    const RT& Bnext_x, const RT& Bnext_y, const RT& Bnext_z)
// {
//     const Vec3& M = {M_x,M_y,M_z};
//     const Vec3& B = {Bprev_x,Bprev_y,Bprev_z};

//     /* The anhysteretic magnetization */
//     GpuArray3 M_an;
//     /* Jacobian matrix chi_an = dM_an/dH_e. The anhysteretic contribution to chi */
//     GpuArray6 chi_an;

//     Vec3 H_e = amrex::Dim3_minus(amrex::Dim3_mult(PhysConst::mu0inv, B),amrex::Dim3_mult(-m_oneMinusAlpha, M));

//     /* Compute M_an and chi_an */
//     JAModel::ComputeAnhystereticContribution(M_an, chi_an, H_e);

//     /* The change in the effective field */
//     Dim3 dH_e;
//     Dim3 dBovermu0 = amrex::Dim3_mult(mu0inv, dB); // mu0^{-1} dB
//     Dim3 M_diff = amrex::Dim3_minus(M_an, M); // M_an - M

//     amrex::Real M_diff_mag = amrex::Dim3_norm2(M_diff); // |M_an - M|
//     amrex::Real denominator = m_ja_model_parameters.k() * M_diff_mag; // k |M_an - M|

//     if (denominator != 0) {
//         // Jacobian matrix chi = dM/dH_e = chi_an + chi_irr
//         amrex::Array<amrex::Real, 6> chi(chi_an);
//         // Add the irreversible contribution chi_irr.
//         chi[0] += M_diff[0] * M_diff[0] / denominator; // xx
// // #if AMREX_SPACEDIM > 1
//         chi[1] += M_diff[1] * M_diff[0] / denominator; // xy
//         chi[2] += M_diff[1] * M_diff[1] / denominator; // yy
// // #endif
// // #if AMREX_SPACEDIM > 2
//         chi[3] += M_diff[2] * M_diff[0] / denominator; // xz
//         chi[4] += M_diff[2] * M_diff[1] / denominator; // yz
//         chi[5] += M_diff[2] * M_diff[2] / denominator; // zz
// // #endif
//         /* Compute dH_e using the total contribution */
//         Compute_dHe(dHe, chi, dBovermu0, oneMinusAlpha);
//         /* Check the sign of (M_an - M) dot dH_e */
//         if (M_diff.dotProduct(dH_e) < 0) {
//             /* Use only anhysteretic part of chi if (M_an - M) dot dH_e < 0 */
//             Compute_dHe(dHe, chi_an, dBovermu0, oneMinusAlpha);
//         }
//     }
//     else {
//         Compute_dHe(dHe, chi_an, dBovermu0, oneMinusAlpha);
//     }

//     Dim3 a = {0,0,0};
//     Dim3 b = {1,1,1};

//     Dim3 c = amrex::Dim3_plus(a, b);

//     /* Update the magnetization */
//     M = M + (1._rt / oneMinusAlpha) * (dBovermu0 - dH_e);


//     M = amrex::Dim3_plus(M,amrex::Dim3_mult(1._rt / oneMinusAlpha, amrex::Dim3_minus(dBovermu0, dH_e)));
// }


// void
// JAModel::ComputeAnhystereticContribution (      JAModel::GpuArray3 &M_an,
//                                                 JAModel::GpuArray6 &chi_an,
//                                           const JAModel::GpuArray3 &H_e)
// {
//     using RT = amrex::Real;
//     RT a = m_ja_model_parameters.a();
//     RT c_Ms_over_a = m_ja_model_parameters.c_times_Ms_over_a();
//     RT H_e_mag = H_e.norm(); // |H_e|
//     // if (H_e_mag <= amrex::Real(0.0)) { return; }
//     RT x = H_e_mag / a; // x = |H_e| / a
//     RT L_over_x = Langevin::L_over_x(x); // L(x)/x
//     RT delta = Langevin::dLdx_minus_L_over_x_over_x2(x)/(a*a);  // delta = L'(x) - L(x)/x

//     M_an = (m_ja_model_parameters.Ms_over_a() * L_over_x) * H_e;

//     // chi_an = dM_an / dH_e is a symmetric matrix, so only 6 components need to be specified.
//     // The order of elements is xx, xy, yy, xz, yz, zz)
//     chi_an[0] = c_Ms_over_a * (L_over_x + delta * H_e[0] * H_e[0]);
//     chi_an[1] = c_Ms_over_a * (           delta * H_e[0] * H_e[1]);
//     chi_an[2] = c_Ms_over_a * (L_over_x + delta * H_e[1] * H_e[1]);
//     chi_an[3] = c_Ms_over_a * (           delta * H_e[0] * H_e[2]);
//     chi_an[4] = c_Ms_over_a * (           delta * H_e[1] * H_e[2]);
//     chi_an[5] = c_Ms_over_a * (L_over_x + delta * H_e[2] * H_e[2]);
// }


// void
// JAModel::ComputeAnhystereticContribution (GpuArray3& M_an, GpuArray6& chi_an, const GpuArray3& H_e)
// {
//     amrex::Real a = m_ja_model_parameters.a();
//     /* c Ms / a */
//     amrex::Real c_Ms_over_a = m_ja_model_parameters.c_times_Ms_over_a();
//     amrex::Real H_e_mag = H_e.vectorLength(); // |H_e|
//     if (H_e_mag <= amrex::Real(0.0)) { return; }
//     amrex::Real x = H_e_mag / a; // x = |H_e| / a
//     amrex::Real L_over_x = Langevin::L_over_x(x); // L(x)/x
//     amrex::Real dLdx = Langevin::dLdx(x); // L'(x)
//     M_an = ((ja_model_parameters->Ms_over_a()) * L_over_x) * H_e;
//     amrex::Real delta = dLdx - L_over_x;  // delta = L'(x) - L(x)/x
//     amrex::RealVect H_e_hat = H_e / H_e_mag; // Unit vector in direction of H_e
//     // chi_an = dM_an / dH_e is a symmetric matrix, so only 6 components need to be specified.
//     // The order of elements is xx, xy, yy, xz, yz, zz)
//     chi_an[0] = c_Ms_over_a * (L_over_x + delta * H_e[0] * H_e[0]);
// #if AMREX_SPACEDIM > 1
//     chi_an[1] = c_Ms_over_a * (           delta * H_e_hat[1] * H_e_hat[0]);
//     chi_an[2] = c_Ms_over_a * (L_over_x + delta * H_e_hat[1] * H_e_hat[1]);
// #endif
// #if AMREX_SPACEDIM > 2
//     chi_an[3] = c_Ms_over_a * (           delta * H_e_hat[2] * H_e_hat[0]);
//     chi_an[4] = c_Ms_over_a * (           delta * H_e_hat[2] * H_e_hat[1]);
//     chi_an[5] = c_Ms_over_a * (L_over_x + delta * H_e_hat[2] * H_e_hat[2]);
// #endif
// }

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

