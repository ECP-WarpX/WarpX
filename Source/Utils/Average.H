#ifndef WARPX_AVERAGE_H_
#define WARPX_AVERAGE_H_

#include "WarpX.H"

namespace Average{

    using namespace amrex;

    /**
     * \brief Returns the cell-centered average of the floating point data contained
     *        in the input Array4 \c mf_in_arr.
     *
     * \param[in] mf_in_arr floating point data to be averaged
     * \param[in] stag      staggering (index type) of the data
     * \param[in] i         index along x to access the Array4 \c mf_in_arr containing
     *                      the data to be averaged
     * \param[in] j         index along y to access the Array4 \c mf_in_arr containing
     *                      the data to be averaged
     * \param[in] k         index along z to access the Array4 \c mf_in_arr containing
     *                      the data to be averaged
     * \param[in] comp      index along the fourth component of the Array4 \c mf_in_arr
     *                      containing the data to be averaged
     *
     * \return averaged field at cell (i,j,k)
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Real ToCellCenter ( Array4<Real const> const& mf_in_arr,
                        const IntVect stag,
                        const int i,
                        const int j,
                        const int k,
                        const int comp=0 )
    {
        const int sx = stag[0];
        const int sy = stag[1];
#if   (AMREX_SPACEDIM == 2)
        constexpr int sz = 0;
#elif (AMREX_SPACEDIM == 3)
        const int sz = stag[2];
#endif
        return 0.125_rt * (   mf_in_arr(i   ,j   ,k   ,comp)
                            + mf_in_arr(i+sx,j   ,k   ,comp)
                            + mf_in_arr(i   ,j+sy,k   ,comp)
                            + mf_in_arr(i   ,j   ,k+sz,comp)
                            + mf_in_arr(i+sx,j+sy,k   ,comp)
                            + mf_in_arr(i   ,j+sy,k+sz,comp)
                            + mf_in_arr(i+sx,j   ,k+sz,comp)
                            + mf_in_arr(i+sx,j+sy,k+sz,comp) );
    };

    /**
     * \brief Stores the cell-centered average of the floating point data contained
     *        in the input MultiFab \c mf_in into the output MultiFab \c mf_out.
     *
     * \param[in,out] mf_out MultiFab containing the floating point data to be filled
     *                       with cell-centered averages
     * \param[in]     mf_in  MultiFab containing the floating point data to be averaged
     * \param[in]     dcomp  offset for the fourth component of the Array4 object,
     *                       extracted from its MultiFab, where the cell-centered
     *                       averages will be stored
     * \param[in]     scomp  optional offset for the fourth component of the Array4
     *                       object, extracted from its MultiFab, containing the data
     *                       to be averaged
     * \param[in]     ncomp  optional number of components to loop over for both the
     *                       input and the output Array4 objects, extracted from their
     *                       respective MultiFab
     */
    void ToCellCenter ( MultiFab& mf_out,
                        const MultiFab& mf_in,
                        const int dcomp,
                        const int ngrow,
                        const int scomp=0,
                        const int ncomp=1 );

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Real CoarsenAndInterpolateKernel ( Array4<Real const> const& mf_fp_arr,
                                       const int sf[3],
                                       const int sc[3],
                                       const int cr[3],
                                       const int i,
                                       const int j,
                                       const int k,
                                       const int comp )
    {
        // Indices of coarse array
        const int ic[3] = { i, j, k };

        // Auxiliary integer arrays
        int num_pts[3];
        int idx_min[3];

        // Compute number of points to loop over (either 1 or 2)
        // and starting indices of fine array in each direction
        for ( int l = 0; l < 3; ++l ) {
            if ( cr[l] == 1 ) {
                num_pts[l] = 1+abs(sf[l]-sc[l]);
                idx_min[l] = ic[l]-sc[l]*(1-sf[l]);
            }
            else {
                num_pts[l] = 2-sf[l];
                idx_min[l] = ic[l]*cr[l]+static_cast<int>(cr[l]/2)*(1-sc[l])-(1-sf[l]);
            }
        }

        // Auxiliary integer variables
        const int numx = num_pts[0];
        const int numy = num_pts[1];
        const int numz = num_pts[2];
        const int imin = idx_min[0];
        const int jmin = idx_min[1];
        const int kmin = idx_min[2];

        // Interpolate over points computed above
        Real c = 0.0;
        for (int kref = 0; kref < numz; ++kref) {
            for (int jref = 0; jref < numy; ++jref) {
                for (int iref = 0; iref < numx; ++iref) {
                    c += mf_fp_arr(imin+iref,jmin+jref,kmin+kref,comp);
                }
            }
        }
        return c / (numx*numy*numz);
    };

    void CoarsenAndInterpolateLoop ( MultiFab& mf_cp,
                                     const MultiFab& mf_fp,
                                     const int dcomp,
                                     const int scomp,
                                     const int ncomp,
                                     const IntVect ratio );

    void CoarsenAndInterpolate ( MultiFab& mf_cp,
                                 const MultiFab& mf_fp,
                                 const int dcomp,
                                 const int scomp,
                                 const int ncomp,
                                 const IntVect ratio );
}

#endif // WARPX_AVERAGE_H_
