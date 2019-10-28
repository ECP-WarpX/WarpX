#ifndef WARPX_SPECTRALHANKELTRANSFORMER_H_
#define WARPX_SPECTRALHANKELTRANSFORMER_H_

#include <AMReX_FArrayBox.H>
#include <HankelTransform.H>

/* \brief Object that allows to transform the fields back and forth between the
 *  spectral and interpolation grid.
 *
 *  Attributes :
 *  - dht0, dhtm, dhtp : the discrete Hankel transform objects
 *     that operates along r
*/

class SpectralHankelTransformer
{

    public:

        SpectralHankelTransformer() :n_rz_azimuthal_modes(0) {};

        SpectralHankelTransformer(const int nr,
                                  const int n_rz_azimuthal_modes,
                                  const amrex::Real rmax);

        // Converts a scalar field from the spectral to the spectral space
        void
        PhysicalToSpectral_Scalar(amrex::Box const & box,
                                  amrex::FArrayBox const & F_physical,
                                  amrex::FArrayBox const & G_spectral);

        // Converts a vector field from the spectral to the spectral space
        void
        PhysicalToSpectral_Vector(amrex::Box const & box,
                                  amrex::FArrayBox const & F_r_physical,
                                  amrex::FArrayBox const & F_t_physical,
                                  amrex::FArrayBox const & G_p_spectral,
                                  amrex::FArrayBox const & G_m_spectral);

        // Converts a scalar field from the spectral to the physical space
        void
        SpectralToPhysical_Scalar(amrex::Box const & box,
                                  amrex::FArrayBox const & G_spectral,
                                  amrex::FArrayBox const & F_physical);

        // Converts a vector field from the spectral to the physical space
        void
        SpectralToPhysical_Vector(amrex::Box const & box,
                                  amrex::FArrayBox const & G_p_spectral,
                                  amrex::FArrayBox const & G_m_spectral,
                                  amrex::FArrayBox const & F_r_physical,
                                  amrex::FArrayBox const & F_t_physical);

    private:

        int n_rz_azimuthal_modes;

        amrex::Vector<HankelTransform> dht0;
        amrex::Vector<HankelTransform> dhtm;
        amrex::Vector<HankelTransform> dhtp;

};

#endif
