/* Copyright 2026
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SCATTERING_UTILS_H_
#define WARPX_SCATTERING_UTILS_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace ScatteringUtils {

    struct AnisotropicCoefficientTable
    {
        amrex::ParticleReal const* m_energies = nullptr;     ///< Sorted center-of-mass energy grid [MeV]
        amrex::ParticleReal const* m_coefficients = nullptr; ///< Row-major coefficients: row i holds m_num_coefficients values for energy i
        int m_num_energies = 0;      ///< Number of rows (energy points) in the table
        int m_num_coefficients = 0;  ///< Number of coefficients per energy row

        /**
         * \brief Return a linearly interpolated coefficient for the given energy.
         *
         * When \p energy_MeV is below the minimum or above the maximum table energy,
         * the coefficient at the nearest endpoint row is returned unchanged (constant
         * extrapolation / clamping). This is a deliberate, conservative choice that
         * keeps the angular distribution physical, but it silently discards any
         * energy dependence outside the table range. The caller checks the energy
         * range once per collision event and reports any out-of-range energies after
         * the particle kernel has completed.
         *
         * @param[in] energy_MeV centre-of-mass energy in MeV
         * @param[in] coefficient_index column index into the coefficient table
         * @return interpolated (or clamped) coefficient value
         */
        [[nodiscard]]
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        amrex::ParticleReal getCoefficient (
            amrex::ParticleReal const energy_MeV,
            int const coefficient_index) const
        {
            using namespace amrex::literals;

            if (energy_MeV <= m_energies[0]) {
                return m_coefficients[coefficient_index];
            }
            if (energy_MeV >= m_energies[m_num_energies - 1]) {
                return m_coefficients[(m_num_energies - 1) * m_num_coefficients
                                      + coefficient_index];
            }

            // Linear scan to find the bracketing interval [lo, hi]
            // TODO: replace with a binary search (O(log N)) since m_energies is sorted
            int hi = 1;
            while ((hi < m_num_energies - 1) && (energy_MeV > m_energies[hi])) {
                ++hi;
            }
            int const lo = hi - 1;

            amrex::ParticleReal const energy_lo = m_energies[lo];
            amrex::ParticleReal const energy_hi = m_energies[hi];
            amrex::ParticleReal const coeff_lo = m_coefficients[lo * m_num_coefficients + coefficient_index];
            amrex::ParticleReal const coeff_hi = m_coefficients[hi * m_num_coefficients + coefficient_index];
            amrex::ParticleReal const interp = (energy_MeV - energy_lo) / (energy_hi - energy_lo); // normalised distance in [0,1]

            return coeff_lo + (coeff_hi - coeff_lo) * interp; // linear interpolation
        }
    };

    /**
     * \brief Sample the cosine of the fusion scattering angle from a coefficient
     * table, linearly interpolating the coefficients in energy.
     *
     * The probability distribution is expressed in terms of Legendre coefficients \f$L_l\f$ as
     * \f$f(\mu) = 1/2 + \sum_{l=1}^{l_{\mathrm{max}}} (2l+1) L_l P_l(\mu)/2\f$, with \f$\mu = \cos(\theta)\f$.
     *
     * @param[in] a_Ll table of anisotropic Legendre coefficients
     * @param[in] a_energy_MeV center-of-mass kinetic energy in MeV
     * @param[in] a_U uniform random number in [0,1]
     * @param[out] a_energy_range_status per-event status: -1 below, 0 within, 1 above
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::ParticleReal sampleCosTheta (
        AnisotropicCoefficientTable const& a_Ll,
        amrex::ParticleReal const a_energy_MeV,
        amrex::ParticleReal const a_U,
        int* const a_energy_range_status)
    {
        using namespace amrex::literals;

        // TODO: In single precision the tolerance might need to be adjusted.
        constexpr amrex::ParticleReal tol = 1.0e-8_prt;
        constexpr int iter_max = 100;
        int const l_max = a_Ll.m_num_coefficients - 1;

        // This status is private to the current collision event. Recording it here avoids
        // repeating the range check for every coefficient and every Newton iteration.
        if (a_energy_range_status != nullptr) {
            if (a_energy_MeV < a_Ll.m_energies[0]) {
                *a_energy_range_status = -1;
            } else if (a_energy_MeV > a_Ll.m_energies[a_Ll.m_num_energies - 1]) {
                *a_energy_range_status = 1;
            }
        }

        amrex::ParticleReal mu = 0.0_prt;

        for (int iter = 0; iter < iter_max; ++iter) {
            amrex::ParticleReal fmu = 0.5_prt;
            amrex::ParticleReal Fmu = 0.5_prt * (mu + 1.0_prt);

            // Generate P_l(mu) in place to keep this routine GPU-callable.
            amrex::ParticleReal P_l_minus_one = 1.0_prt;
            amrex::ParticleReal P_l = mu;
            for (int l = 1; l <= l_max; ++l) {
                auto const l_prt = static_cast<amrex::ParticleReal>(l);
                amrex::ParticleReal const P_l_plus_one =
                    ((2.0_prt * l_prt + 1.0_prt) * mu * P_l
                     - l_prt * P_l_minus_one) / (l_prt + 1.0_prt);
                amrex::ParticleReal const coefficient =
                    a_Ll.getCoefficient(a_energy_MeV, l);

                fmu += (l_prt + 0.5_prt) * coefficient * P_l;
                Fmu += 0.5_prt * coefficient * (P_l_plus_one - P_l_minus_one);

                P_l_minus_one = P_l;
                P_l = P_l_plus_one;
            }

            amrex::ParticleReal const residual = Fmu - a_U;
            if (std::abs(residual) < tol) {
                return mu;
            }

            amrex::ParticleReal const dmu = residual / fmu;
            mu -= dmu;
            mu = amrex::max(-1.0_prt, amrex::min(1.0_prt, mu));

            if (std::abs(dmu) < tol) {
                return mu;
            }
        }

        return mu;
    }
}

#endif // WARPX_SCATTERING_UTILS_H_
