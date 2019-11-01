#ifndef WARPX_SPECTRAL_HANKEL_K_SPACE_H_
#define WARPX_SPECTRAL_HANKEL_K_SPACE_H_

#include <SpectralKSpace.H>

/* \brief Class that represents the spectral, Hankel/FFT, space.
 *
 * (Contains info about the size of the spectral space corresponding
 * to each box in `realspace_ba`, as well as the value of the
 * corresponding k coordinates)
 */
class SpectralHankelKSpace
    :
    public SpectralKSpace
{
    public:
        SpectralHankelKSpace(const amrex::BoxArray& realspace_ba,
                             const amrex::DistributionMapping& dm,
                             const amrex::RealVect realspace_dx);
};

#endif // WARPX_SPECTRAL_HANKEL_K_SPACE_H_
