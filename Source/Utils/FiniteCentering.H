#ifndef WARPX_FINITE_CENTERING_H_
#define WARPX_FINITE_CENTERING_H_

#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

namespace FiniteCentering{

    using namespace amrex;

    // TODO Add Doxygen
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE

    void compute_centering_coefficients(amrex::GpuArray<amrex::Real, 8>& cc,
                                        const amrex::Real r, const int m)
    {
        if (m == 4)
        {
            const amrex::Real r2 = r*r;
            const amrex::Real r3 = r*r2;
            cc[0] = - r3/6. + r/6.;
            cc[1] = r3/2. + r2/2. - r;
            cc[2] = - r3/2. - r2 + r/2. + 1.;
            cc[3] = r3/6. + r2/2. + r/3.;
        }
        else if (m == 6)
        {
            const amrex::Real r2 = r*r;
            const amrex::Real r3 = r*r2;
            const amrex::Real r4 = r*r3;
            const amrex::Real r5 = r*r4;
            cc[0] = - r5/120. + r3/24. - r/30.;
            cc[1] = r5/24. + r4/24. - r3*7./24. - r2/24. + r/4.;
            cc[2] = - r5/12. - r4/6. + r3*7./12. + r2*2./3. - r;
            cc[3] = r5/12. + r4/4. - r3*5./12. - r2*5./4. + r/3. + 1.;
            cc[4] = - r5/24. - r4/6. + r3/24. + r2*2./3. + r/2.;
            cc[5] = r5/120. + r4/24. + r3/24. - r2/24. - r/20.;
        }
        else if (m == 8)
        {
            const amrex::Real r2 = r*r;
            const amrex::Real r3 = r*r2;
            const amrex::Real r4 = r*r3;
            const amrex::Real r5 = r*r4;
            const amrex::Real r6 = r*r5;
            const amrex::Real r7 = r*r6;
            cc[0] = - r7/5040. + r5/360. - r3*7./720. + r/140.;
            cc[1] = r7/720. + r6/720. - r5*17./720. - r4/144. + r3*4./45. + r2/180. - r/15.;
            cc[2] = - r7/240. - r6/120. + r5*3./40. + r4/12. - r3*89./240. - r2*3./40. + r*3./10.;
            cc[3] = r7/144. + r6/48. - r5*17./144. - r4*13./48. + r3*11./18. + r2*3./4. - r;
            cc[4] = - r7/144. - r6/36. + r5*7./72. + r4*7./18. - r3*49./144. - r2*49./36. + r/4. + 1.;
            cc[5] = r7/240. + r6/48. - r5*3./80. - r4*13./48. - r3/15. + r2*3./4. + r*3./5.;
            cc[6] = - r7/720. - r6/120. + r5/360. + r4/12. + r3*71./720. - r2*3./40. - r/10.;
            cc[7] = r7/5040. + r6/720. + r5/720. - r4/144. - r3/90. + r2/180. + r/105.;
        }
        else
        {
            amrex::Abort("Centering coefficients implemented only for orders 4,6,8");
        }
    }
}

#endif // WARPX_FINITE_CENTERING_H_
