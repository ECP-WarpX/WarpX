#include "Interpolate.H"
#include "Interpolate_K.H"

namespace Interpolate
{
    using namespace amrex;

    std::unique_ptr<MultiFab>
    getInterpolatedScalar(
        const MultiFab& F_cp, const MultiFab& F_fp,
        const DistributionMapping& dm, const int r_ratio,
        const Real* /*dx*/, const int ngrow )
    {
        // Prepare the structure that will contain the returned fields
        std::unique_ptr<MultiFab> interpolated_F;
        interpolated_F.reset( new MultiFab(F_fp.boxArray(), dm, 1, ngrow) );
        interpolated_F->setVal(0.);

        // Loop through the boxes and interpolate the values from the _cp data
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox ffab; // Temporary array ; contains interpolated fields
            for (MFIter mfi(*interpolated_F); mfi.isValid(); ++mfi)
            {
                Box finebx = mfi.fabbox();
                finebx.coarsen(r_ratio).refine(r_ratio); // so that finebx is coarsenable

                const FArrayBox& cfab = (F_cp)[mfi];
                ffab.resize(finebx);

                // - Fully nodal
                if ( F_fp.is_nodal() ){
                    IntVect refinement_vector{AMREX_D_DECL(r_ratio, r_ratio, r_ratio)};
                    node_bilinear_interp.interp(cfab, 0, ffab, 0, 1,
                                                finebx, refinement_vector, {}, {}, {}, 0, 0, RunOn::Cpu);
                } else {
                    amrex::Abort("Unknown field staggering.");
                }

                // Add temporary array to the returned structure
                const Box& bx = (*interpolated_F)[mfi].box();
                (*interpolated_F)[mfi].plus<RunOn::Host>(ffab, bx, bx, 0, 0, 1);
            }
        }
        return interpolated_F;
    }

    std::array<std::unique_ptr<MultiFab>, 3>
    getInterpolatedVector(
        const MultiFab* Fx_cp,
        const MultiFab* Fy_cp,
        const MultiFab* Fz_cp,
        const MultiFab* Fx_fp,
        const MultiFab* Fy_fp,
        const MultiFab* Fz_fp,
        const DistributionMapping& dm, const int r_ratio,
        const Real* dx, const int ngrow )
    {

        // Prepare the structure that will contain the returned fields
        std::array<std::unique_ptr<MultiFab>, 3> interpolated_F;
        interpolated_F[0].reset( new MultiFab(Fx_fp->boxArray(), dm, 1, ngrow) );
        interpolated_F[1].reset( new MultiFab(Fy_fp->boxArray(), dm, 1, ngrow) );
        interpolated_F[2].reset( new MultiFab(Fz_fp->boxArray(), dm, 1, ngrow) );

        IntVect fx_type = interpolated_F[0]->ixType().toIntVect();
        IntVect fy_type = interpolated_F[1]->ixType().toIntVect();
        IntVect fz_type = interpolated_F[2]->ixType().toIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*interpolated_F[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& boxx = mfi.growntilebox(fx_type);
            Box const& boxy = mfi.growntilebox(fy_type);
            Box const& boxz = mfi.growntilebox(fz_type);

            Array4<Real      > const& fx = interpolated_F[0]->array(mfi);
            Array4<Real      > const& fy = interpolated_F[1]->array(mfi);
            Array4<Real      > const& fz = interpolated_F[2]->array(mfi);
            Array4<Real const> const& cx = Fx_cp->const_array(mfi);
            Array4<Real const> const& cy = Fy_cp->const_array(mfi);
            Array4<Real const> const& cz = Fz_cp->const_array(mfi);

            amrex::ParallelFor(boxx, boxy, boxz,
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    interp(j,k,l,fx,cx,r_ratio,fx_type);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    interp(j,k,l,fy,cy,r_ratio,fy_type);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    interp(j,k,l,fz,cz,r_ratio,fz_type);
                });
        }
        return interpolated_F;
    }
}
