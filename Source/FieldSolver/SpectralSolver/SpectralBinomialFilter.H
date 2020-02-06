#ifndef WARPX_SPECTRAL_BINOMIAL_FILTER_H_
#define WARPX_SPECTRAL_BINOMIAL_FILTER_H_

#include <HankelTransform.H>
#include <SpectralKSpace.H>

using namespace amrex;

// `KFilterArray` holds a one 1D array ("ManagedVector") that
// implements the filter.
using KFilterArray = Gpu::ManagedVector<Real>;

/**
 * \brief Class that sets up binomial filtering in k space
 *
 * Contains the filter arrays for r and z.
 */
class SpectralBinomialFilter
{
    public:
        SpectralBinomialFilter() {};
        void InitFilterArrayR(HankelTransform::RealVector const & kr,
                              Real const dr,
                              int const npasses,
                              bool const compensation);
        void InitFilterArrayZ(RealKVector const & kz,
                              Real const dz,
                              int const npasses,
                              bool const compensation);

        KFilterArray const & getFilterArrayR() {return filter_r;};
        KFilterArray const & getFilterArrayZ() {return filter_z;};

    protected:

        KFilterArray filter_r;
        KFilterArray filter_z;

};

#endif
