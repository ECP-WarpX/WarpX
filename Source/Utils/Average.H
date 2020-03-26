#ifndef WARPX_AVERAGE_H_
#define WARPX_AVERAGE_H_

#include "WarpX.H"

namespace Average{

    using namespace amrex;

    /**
     * \brief Returns the cell-centered average of the floating point data contained
     *        in the input
     *        <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *        \c mf_in_arr.
     *
     * \param[in] mf_in_arr floating point data to be averaged
     * \param[in] stag      staggering (index type) of the data
     * \param[in] i         index along x to access the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
     * \param[in] j         index along y to access the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
     * \param[in] k         index along z to access the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
     * \param[in] comp      index along the fourth component of the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
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
     *        in the input
     *        <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *        \c mf_in into the output
     *        <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *        \c mf_out.
     *
     * \param[in,out] mf_out <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       containing the floating point data to be filled with cell-centered averages
     * \param[in]     mf_in  <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       containing the floating point data to be averaged
     * \param[in]     dcomp  offset for the fourth component of the
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                       object, extracted from its
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       , where the cell-centered averages will be stored
     * \param[in]     scomp  optional offset for the fourth component of the
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                       object, extracted from its
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       , containing the data to be averaged
     * \param[in]     ncomp  optional number of components to loop over for both the
     *                       input and the output
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                       objects, extracted from their respective
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     */
    void ToCellCenter ( MultiFab& mf_out,
                        const MultiFab& mf_in,
                        const int dcomp,
                        const int ngrow,
                        const int scomp=0,
                        const int ncomp=1 );

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    void FineToCoarse ( Box const& bx,
                        Array4<Real> const& crse,
                        Array4<Real const> const& fine,
                        const IntVect stag,
                        int ccomp,
                        int fcomp,
                        int ncomp,
                        IntVect const& ratio )
    {

        const int sx = stag[0];
        const int sy = stag[1];
#if   (AMREX_SPACEDIM == 2)
        constexpr int sz = 0;
#elif (AMREX_SPACEDIM == 3)
        const int sz = stag[2];
#endif
        const auto clo = lbound(bx);
        const auto chi = ubound(bx);
        const int facx = ratio[0];
        const int facy = ratio[1];
#if   (AMREX_SPACEDIM == 2)
        const int facz = 1;
#elif (AMREX_SPACEDIM == 3)
        const int facz = ratio[2];
#endif
        const int itx = facx*(1-sx)+sx;
        const int ity = facy*(1-sy)+sy;
        const int itz = facz*(1-sz)+sz;
        const Real volfrac = 1.0/static_cast<Real>(itx*ity*itz);

        for (int n = 0; n < ncomp; ++n) {
            for (int k = clo.z; k <= chi.z; ++k) {
                int kk = k*facz;
                for (int j = clo.y; j <= chi.y; ++j) {
                    int jj = j*facy;
                    for (int i = clo.x; i <= chi.x; ++i) {
                        int ii = i*facx;
                        Real c = 0.;
                        for (int kref = 0; kref < itz; ++kref) {
                            for (int jref = 0; jref < ity; ++jref) {
                                for (int iref = 0; iref < itx; ++iref) {
                                    c += fine(ii+iref,jj+jref,kk+kref,n+fcomp);
                                }
                            }
                        }
                        crse(i,j,k,n+ccomp) = volfrac * c;
                    }
                }
            }
        }
    };

    void FineToCoarse ( const MultiFab& S_fine,
                        MultiFab& S_crse,
                        int scomp,
                        int ncomp,
                        const IntVect& ratio );
}

#endif // WARPX_AVERAGE_H_
