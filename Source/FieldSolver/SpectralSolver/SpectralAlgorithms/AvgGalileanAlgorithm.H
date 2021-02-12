#ifndef WARPX_AVG_GALILEAN_ALGORITHM_H_
#define WARPX_AVG_GALILEAN_ALGORITHM_H_

#include "FieldSolver/SpectralSolver/SpectralAlgorithms/SpectralBaseAlgorithm.H"

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class AvgGalileanAlgorithm : public SpectralBaseAlgorithm
{
    public:
        AvgGalileanAlgorithm (const SpectralKSpace& spectral_kspace,
                              const amrex::DistributionMapping& dm,
                              const int norder_x, const int norder_y,
                              const int norder_z, const bool nodal,
                              const amrex::Array<amrex::Real,3>& v_galilean,
                              const amrex::Real dt);
        // Redefine update equation from base class
        virtual void pushSpectralFields (SpectralFieldData& f) const override final;
        virtual int getRequiredNumberOfFields () const override final {
             return SpectralAvgFieldIndex::n_fields;
        }
        void InitializeSpectralCoefficients(
           const SpectralKSpace& spectral_kspace,
           const amrex::DistributionMapping& dm,
           const amrex::Array<amrex::Real, 3>& v_galilean,
           const amrex::Real dt);

        /**
         * \brief Virtual function for current correction in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c CurrentCorrection in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         * \param[in]     rho        Unique pointer to \c MultiFab storing the charge density
         */
        virtual void CurrentCorrection (SpectralFieldData& field_data,
                                        std::array<std::unique_ptr<amrex::MultiFab>,3>& current,
                                        const std::unique_ptr<amrex::MultiFab>& rho) override final;

        /**
         * \brief Virtual function for Vay current deposition in Fourier space
         * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>).
         * This function overrides the virtual function \c VayDeposition in the
         * base class \c SpectralBaseAlgorithm and cannot be overridden by further
         * derived classes.
         *
         * \param[in,out] field_data All fields in Fourier space
         * \param[in,out] current    Array of unique pointers to \c MultiFab storing
         *                           the three components of the current density
         */
        virtual void VayDeposition (SpectralFieldData& field_data,
                                    std::array<std::unique_ptr<amrex::MultiFab>,3>& current) override final;

    private:
        SpectralRealCoefficients C_coef, S_ck_coef, C1_coef, C3_coef, S1_coef,S3_coef;
        SpectralComplexCoefficients Theta2_coef, X1_coef, X2_coef, X3_coef, X4_coef, Psi1_coef, Psi2_coef, Psi3_coef, Psi4_coef, A1_coef, A2_coef, Rhoold_coef, Rhonew_coef, Jcoef_coef;

};

#endif // WARPX_GALILEAN_ALGORITHM_H_
