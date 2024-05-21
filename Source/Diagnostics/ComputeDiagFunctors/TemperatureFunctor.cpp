
#include "TemperatureFunctor.H"

#include "Diagnostics/ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>

#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

TemperatureFunctor::TemperatureFunctor (const int lev,
        const amrex::IntVect crse_ratio, const int ispec, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_lev(lev), m_ispec(ispec)
{
    // Write only in one output component.
    AMREX_ALWAYS_ASSERT(ncomp == 1);
}

void
TemperatureFunctor::operator() (amrex::MultiFab& mf_dst, const int dcomp, const int /*i_buffer*/) const
{
    using namespace amrex::literals;
    auto& warpx = WarpX::GetInstance();

    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;

    // Temporary cell-centered, multi-component MultiFab for storing particles sums and result
    amrex::MultiFab sum_mf(warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 7, ng);

    auto& pc = warpx.GetPartContainer().GetParticleContainer(m_ispec);
    amrex::Real mass = pc.getMass();  // Note, implicit conversion from ParticleReal

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(mass > 0.,
        "The temperature diagnostic can not be calculated for a massless species.");

    ParticleToMesh(pc, sum_mf, m_lev,
            [=] AMREX_GPU_DEVICE (const WarpXParticleContainer::SuperParticleType& p,
                amrex::Array4<amrex::Real> const& out_array,
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
            {
                // Get position in AMReX convention to calculate corresponding index.
                // Ideally this will be replaced with the AMReX NGP interpolator
                // Always do x direction.
                int ii = 0, jj = 0, kk = 0;
                const amrex::ParticleReal x = p.pos(0);
                const amrex::Real lx = (x - plo[0]) * dxi[0];
                ii = static_cast<int>(amrex::Math::floor(lx));
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
                const amrex::ParticleReal y = p.pos(1);
                const amrex::Real ly = (y - plo[1]) * dxi[1];
                jj = static_cast<int>(amrex::Math::floor(ly));
#endif
#if defined(WARPX_DIM_3D)
                const amrex::ParticleReal z = p.pos(2);
                const amrex::Real lz = (z - plo[2]) * dxi[2];
                kk = static_cast<int>(amrex::Math::floor(lz));
#endif

                const amrex::ParticleReal w  = p.rdata(PIdx::w);
                const amrex::ParticleReal ux = p.rdata(PIdx::ux);
                const amrex::ParticleReal uy = p.rdata(PIdx::uy);
                const amrex::ParticleReal uz = p.rdata(PIdx::uz);
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 0), (amrex::Real)(w));
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 1), (amrex::Real)(w*ux));
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 2), (amrex::Real)(w*uy));
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 3), (amrex::Real)(w*uz));
            });

    // Divide value by number of particles for average
    for (amrex::MFIter mfi(sum_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        amrex::Array4<amrex::Real> const& out_array = sum_mf.array(mfi);
        amrex::ParallelFor(box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (out_array(i,j,k,0) > 0) {
                        const amrex::Real invsum = 1._rt/out_array(i,j,k,0);
                        out_array(i,j,k,1) *= invsum;
                        out_array(i,j,k,2) *= invsum;
                        out_array(i,j,k,3) *= invsum;
                    }
                });
    }

    // Calculate the sum of the squares, subtracting the averages
    // These loops must be written out since ParticleToMesh always zeros out the mf.
    const auto plo = pc.Geom(m_lev).ProbLoArray();
    const auto dxi = pc.Geom(m_lev).InvCellSizeArray();
    for (WarpXParIter pti(pc, m_lev); pti.isValid(); ++pti)
    {
        const long np = pti.numParticles();
        auto const & wp = pti.GetAttribs(PIdx::w);
        auto const & uxp = pti.GetAttribs(PIdx::ux);
        auto const & uyp = pti.GetAttribs(PIdx::uy);
        auto const & uzp = pti.GetAttribs(PIdx::uz);

        auto const GetPosition = GetParticlePosition<PIdx>(pti);

        amrex::Array4<amrex::Real> out_array = sum_mf[pti].array();

        amrex::ParallelFor(np,
            [=] AMREX_GPU_DEVICE (long ip) {
                // --- Get particle quantities
                amrex::ParticleReal xp, yp, zp;
                GetPosition.AsStored(ip, xp, yp, zp);

                // Get position in AMReX convention to calculate corresponding index.
                int ii = 0, jj = 0, kk = 0;
                const amrex::Real lx = (xp - plo[0]) * dxi[0];
                ii = static_cast<int>(amrex::Math::floor(lx));
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const amrex::Real lz = (zp - plo[1]) * dxi[1];
                jj = static_cast<int>(amrex::Math::floor(lz));
#elif defined(WARPX_DIM_3D)
                const amrex::Real ly = (yp - plo[1]) * dxi[1];
                jj = static_cast<int>(amrex::Math::floor(ly));
                const amrex::Real lz = (zp - plo[2]) * dxi[2];
                kk = static_cast<int>(amrex::Math::floor(lz));
#endif

                const amrex::ParticleReal w  = wp[ip];
                const amrex::ParticleReal ux = uxp[ip] - out_array(ii, jj, kk, 1);
                const amrex::ParticleReal uy = uyp[ip] - out_array(ii, jj, kk, 2);
                const amrex::ParticleReal uz = uzp[ip] - out_array(ii, jj, kk, 3);
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 4), (amrex::Real)(w*ux*ux));
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 5), (amrex::Real)(w*uy*uy));
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 6), (amrex::Real)(w*uz*uz));
            });
    }

    // Divide the squares by number of particles for average and calculate the temperature
    for (amrex::MFIter mfi(sum_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        amrex::Array4<amrex::Real> const& out_array = sum_mf.array(mfi);
        amrex::ParallelFor(box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (out_array(i,j,k,0) > 0) {
                        const amrex::Real invsum = 1._rt/out_array(i,j,k,0);
                        out_array(i,j,k,4) *= invsum;
                        out_array(i,j,k,5) *= invsum;
                        out_array(i,j,k,6) *= invsum;
                        out_array(i,j,k,0) = mass*(out_array(i,j,k,4) + out_array(i,j,k,5) + out_array(i,j,k,6))/(3._rt*PhysConst::q_e);
                    }
                });
    }

    // Coarsen and interpolate from sum_mf to the output diagnostic MultiFab, mf_dst.
    ablastr::coarsen::sample::Coarsen(mf_dst, sum_mf, dcomp, 0, nComp(), 0, m_crse_ratio);

}
