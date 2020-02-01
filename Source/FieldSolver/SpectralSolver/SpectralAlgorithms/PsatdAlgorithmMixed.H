#ifndef WARPX_PSATD_ALGORITHM_MIXED_H_
#define WARPX_PSATD_ALGORITHM_MIXED_H_

#include <SpectralBaseAlgorithm.H>

// \brief Class that updates the fields in spectral space and stores the coefficients
// of the corresponding update equations. Objects of this type are instantiated when
// psatd_push_algo = mixed. Implementation of update equation for E (equation 13 of
// reference below) differs from default implementation (class PsatdAlgorithm) as
// only Gauss law is used to express E in terms of rho, while all other terms are
// expressed in terms of J (no use of discrete continuity equation in Fourier space).
// Reference: Vay et al, JCP 243, 2013
class PsatdAlgorithmMixed : public SpectralBaseAlgorithm
{

    public:

        PsatdAlgorithmMixed ( const SpectralKSpace& spectral_kspace,
                              const amrex::DistributionMapping& dm,
                              const int norder_x,
                              const int norder_y,
                              const int norder_z,
                              const bool nodal,
                              const amrex::Real dt );

        // Redefine functions from base class

        virtual void pushSpectralFields ( SpectralFieldData& f ) const override final;

        virtual int getRequiredNumberOfFields () const override final {
            return SpectralFieldIndex::n_fields;
        }

        void InitializeSpectralCoefficients ( const SpectralKSpace& spectral_kspace,
                                              const amrex::DistributionMapping& dm,
                                              const amrex::Real dt );

    private:
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_ALGORITHM_MIXED_H_
