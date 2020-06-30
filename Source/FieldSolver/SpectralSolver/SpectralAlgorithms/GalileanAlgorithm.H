#ifndef WARPX_GALILEAN_ALGORITHM_H_
#define WARPX_GALILEAN_ALGORITHM_H_

#include "SpectralBaseAlgorithm.H"

#if WARPX_USE_PSATD
/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class GalileanAlgorithm : public SpectralBaseAlgorithm
{
    public:
        GalileanAlgorithm(const SpectralKSpace& spectral_kspace,
                         const amrex::DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal,
                         const amrex::Array<amrex::Real,3>& v_galilean,
                         const amrex::Real dt);
        // Redefine update equation from base class
        virtual void pushSpectralFields(SpectralFieldData& f) const override final;
        virtual int getRequiredNumberOfFields() const override final {
            return SpectralFieldIndex::n_fields;
        };
        void InitializeSpectralCoefficients(const SpectralKSpace& spectral_kspace,
                                    const amrex::DistributionMapping& dm,
                                    const amrex::Array<amrex::Real, 3>& v_galilean,
                                    const amrex::Real dt);

        /**
         * \brief Virtual function for current correction in Fourier space
         * (equation (10) of https://doi.org/10.1103/PhysRevE.94.053305).
         * This function overrides the virtual function \c CurrentCorrection in the
         * base class \c SpectralBaseAlgorithm (qualifier \c override) and cannot be
         * overridden by further derived classes (qualifier \c final).
         *
         * \param[in,out] field_data all fields in Fourier space
         * \param[in,out] current    three-dimensional array of unique pointers to MultiFab
         *                           storing the three components of the current density
         * \param[in]     rho        unique pointer to MultiFab storing the charge density
         */
        virtual void CurrentCorrection ( SpectralFieldData& field_data,
                                         std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                         const std::unique_ptr<amrex::MultiFab>& rho ) override final;

    private:
        SpectralRealCoefficients C_coef, S_ck_coef;
        SpectralComplexCoefficients Theta2_coef, X1_coef, X2_coef, X3_coef, X4_coef;
        amrex::Real m_dt;
        amrex::Array<amrex::Real,3> m_v_galilean;
};
#endif // WARPX_USE_PSATD
#endif // WARPX_GALILEAN_ALGORITHM_H_
