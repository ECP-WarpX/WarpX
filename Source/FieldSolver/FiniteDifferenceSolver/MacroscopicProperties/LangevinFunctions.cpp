#include "LangevinFunctions.H"

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