
#include <AMReX_REAL.H>

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