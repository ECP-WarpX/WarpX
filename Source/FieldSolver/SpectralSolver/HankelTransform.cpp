#include <WarpX.H>

#include <HankelTransform.H>
#include <BesselRoots.H>

#include <blas.hh>
#include <lapack.hh>

HankelTransform::HankelTransform (const int hankel_order,
                                  const int azimuthal_mode,
                                  const int nr,
                                  const amrex::Real rmax)
: nr(nr), nk(nr)
{
    // Check that azimuthal_mode has a valid value
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hankel_order-1 <= azimuthal_mode && azimuthal_mode <= hankel_order+1,
                                     "azimuthal_mode must be either hankel_order-1, hankel_order or hankel_order+1");

    RealVector alphas;
    amrex::Vector<int> alpha_errors;

    GetBesselRoots(azimuthal_mode, nk, alphas, alpha_errors);
    // Add of check of alpha_errors, should be all zeros

    // Calculate the spectral grid
    nu.resize(nk);
    for (int ik=0 ; ik < nk ; ik++) {
        nu[ik] = alphas[ik]/(2*MathConst::pi*rmax);
    }

    // Calculate the spatial grid (Uniform grid with a half-cell offset)
    r.resize(nr);
    for (int ir=0 ; ir < nr ; ir++) {
        r[ir] = (rmax/nr)*(ir + 0.5);
    }

    // Calculate and store the inverse matrix invM
    // (imposed by the constraints on the DHT of Bessel modes)
    // NB: When compared with the FBPIC article, all the matrices here
    // are calculated in transposed form. This is done so as to use the
    // `dot` and `gemm` functions, in the `transform` method.
    int p_denom;
    if (hankel_order == azimuthal_mode) {
        p_denom = hankel_order + 1;
    } else {
        p_denom = hankel_order;
    }

    RealVector denom(nk);
    for (int ik=0 ; ik < nk ; ik++) {
        const amrex::Real jna = jn(p_denom, alphas[ik]);
        denom[ik] = MathConst::pi*rmax*rmax*jna*jna;
    }

    RealVector num(nk*nr);
    for (int ir=0 ; ir < nr ; ir++) {
        for (int ik=0 ; ik < nk ; ik++) {
            const int ii = ik + ir*nk;
            num[ii] = jn(hankel_order, 2*MathConst::pi*r[ir]*nu[ik]);
        }
    }

    // Get the inverse matrix
    invM.resize(nk*nr);
    if (azimuthal_mode > 0) {
        for (int ir=0 ; ir < nr ; ir++) {
            for (int ik=1 ; ik < nk ; ik++) {
                const int ii = ik + ir*nk;
                invM[ii] = num[ii]/denom[ik];
            }
        }
        // ik = 0
        // In this case, the functions are represented by Bessel functions
        // *and* an additional mode (below) which satisfies the same
        // algebric relations for curl/div/grad as the regular Bessel modes,
        // with the value kperp=0.
        // The normalization of this mode is arbitrary, and is chosen
        // so that the condition number of invM is close to 1
        if (hankel_order == azimuthal_mode-1) {
            for (int ir=0 ; ir < nr ; ir++) {
                const int ii = ir*nk;
                invM[ii] = std::pow(r[ir], (azimuthal_mode-1))/(MathConst::pi*std::pow(rmax, (azimuthal_mode+1)));
            }
        } else {
            for (int ir=0 ; ir < nr ; ir++) {
                const int ii = ir*nk;
                invM[ii] = 0.;
            }
        }
    } else {
        for (int ir=0 ; ir < nr ; ir++) {
            for (int ik=0 ; ik < nk ; ik++) {
                const int ii = ik + ir*nk;
                invM[ii] = num[ii]/denom[ik];
            }
        }
    }

    // Calculate the matrix M by inverting invM
    if (azimuthal_mode !=0 && hankel_order != azimuthal_mode-1) {

        M.resize(nk*nr, 0.);
        RealVector invMcopy(invM);
        RealVector sdiag(nk-1, 0.);
        RealVector u((nk-1)*(nk-1), 0.);
        RealVector vt((nr)*(nr), 0.);
        RealVector sp((nr)*(nk-1), 0.);
        RealVector temp((nr)*(nk-1), 0.);

        // Note that invMcopy.dataPtr()+1 is passed in so that the first ik row is skipped
        // A copy is passed in since the matrix is destroyed
        lapack::gesvd(lapack::Job::AllVec, lapack::Job::AllVec,
                      nk-1, nr, invMcopy.dataPtr()+1, nk,
                      sdiag.dataPtr(), u.dataPtr(), nk-1,
                      vt.dataPtr(), nr);

        for (int i=0 ; i < nk-1 ; i++) {
            if (sdiag[i] != 0.) {
                const int j = i + i*nk;
                sp[j] = 1./sdiag[i];
            }
        }

        // a_pseudo(1:n,1:m) = matmul(transpose(vt(1:n,1:n)), matmul(sp(1:n,1:m), transpose(u(1:m,1:m))))
        blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans,
                   nr, nk-1, nk-1, 1.,
                   sp.dataPtr(), nr,
                   u.dataPtr(), nk-1, 0.,
                   temp.dataPtr(), nr);
        // Note that M.dataPtr()+nr is passed in so that the first ir column is skipped
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                   nr, nk-1, nr, 1.,
                   vt.dataPtr(), nr,
                   temp.dataPtr(), nr, 0.,
                   M.dataPtr()+nr, nr);

    } else {

        M = invM;
        amrex::Vector<int64_t> ipiv(nr);
        lapack::getrf(nk, nr, M.dataPtr(), nk, ipiv.dataPtr());
        lapack::getri(nr, M.dataPtr(), nr, ipiv.dataPtr());

    }

}

const HankelTransform::RealVector
HankelTransform::getRadialGrid ()
{
    return r;
}


const HankelTransform::RealVector
HankelTransform::getSpectralFrequencies ()
{
    return nu;
}


void
HankelTransform::HankelForwardTransform (const int nz, amrex::Array4<amrex::Real const> const& F, amrex::Array4<amrex::Real> const& G)
{
    // Note that M is flagged to be transposed since it has dimensions (nr, nk)
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
               nk, nz, nr, 1.,
               M.dataPtr(), nk,
               F.dataPtr(), nr, 0.,
               G.dataPtr(), nk);
}


void
HankelTransform::HankelInverseTransform (const int nz, amrex::Array4<amrex::Real const> const& G, amrex::Array4<amrex::Real> const& F)
{
    // Note that invM is flagged to be transposed since it has dimensions (nk, nr)
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
               nr, nz, nk, 1.,
               invM.dataPtr(), nr,
               G.dataPtr(), nk, 0.,
               F.dataPtr(), nr);
}
