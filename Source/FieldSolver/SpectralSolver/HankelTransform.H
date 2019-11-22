#ifndef WARPX_HANKEL_TRANSFORM_H_
#define WARPX_HANKEL_TRANSFORM_H_

/* #include <AMReX_REAL.H> */
/* #include <AMReX_Vector.H> */
/* #include <AMReX_CudaContainers.H> */

/* \brief This defines the class that performs the Hankel transform.
 * Original authors: Remi Lehe, Manuel Kirchen
 *
 * 
 * Definition of the Hankel forward and backward transform of order p:
 * g(\nu) = 2 \pi \int_0^\infty f(r) J_p( 2 \pi \nu r) r dr
 * f( r ) = 2 \pi \int_0^\infty g(\nu) J_p( 2 \pi \nu r) \nu d\nu d
*/
class HankelTransform
{
    public:

        using RealVector = amrex::Gpu::ManagedDeviceVector<amrex::Real>;

        // Constructor
        HankelTransform(const int hankel_order,
                        const int azimuthal_mode,
                        const int nr,
                        const amrex::Real rmax);

        const RealVector & getRadialGrid() {return r;}
        const RealVector & getSpectralFrequencies() {return nu;}

        void HankelForwardTransform(const int nz,
                                    amrex::FArrayBox const& F, int const F_icomp,
                                    amrex::FArrayBox      & G, int const G_icomp);

        void HankelInverseTransform(const int nz,
                                    amrex::FArrayBox const& G, int const G_icomp,
                                    amrex::FArrayBox      & F, int const F_icomp);

    private:
        // Even though nk == nr always, use a seperate variable for clarity.
        int nr, nk;

        RealVector nu;
        RealVector r;

        RealVector invM;
        RealVector M;
};

#endif
