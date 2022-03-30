
#include "ParticleReductionFunctor.H"

#include "Diagnostics/ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/CoarsenIO.H"
#include "WarpX.H"

#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

ParticleReductionFunctor::ParticleReductionFunctor (const amrex::MultiFab* mf_src, const int lev,
        const amrex::IntVect crse_ratio, const std::string fn_str, const int ispec, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_lev(lev), m_ispec(ispec)
{
    // mf_src will not be used, let's make sure it's null.
    AMREX_ALWAYS_ASSERT(mf_src == nullptr);
    // Write only in one output component.
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    // Allocate and compile a parser based on the input string fn_str
    m_map_fn_parser = std::make_unique<amrex::Parser>(makeParser(
                fn_str, {"x", "y", "z", "ux", "uy", "uz"}));
    m_map_fn = m_map_fn_parser->compile<6>();
}

void
ParticleReductionFunctor::operator() (amrex::MultiFab& mf_dst, const int dcomp, const int /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // Temporary cell-centered, multi-component MultiFab for storing particles per cell.
    amrex::MultiFab red_mf(warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 1, ng);
    amrex::MultiFab ppc_mf(warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 1, ng);
    // Set value to 0, and increment the value in each cell with ppc.
    red_mf.setVal(0._rt);
    ppc_mf.setVal(0._rt);
    auto& pc = warpx.GetPartContainer().GetParticleContainer(m_ispec);
    // Copy over compiled parser function so that it can be captured by value in the lambda
    auto map_fn = m_map_fn;
    ParticleToMesh(pc, red_mf, m_lev,
            [=] AMREX_GPU_DEVICE (const WarpXParticleContainer::SuperParticleType& p,
                amrex::Array4<amrex::Real> const& out_array,
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
            {
                // Get position in WarpX convention to use in parser. Will be different from
                // p.pos() for 1D and 2D simulations.
                amrex::ParticleReal xw = 0._rt, yw = 0._rt, zw = 0._rt;
                get_particle_position(p, xw, yw, zw);

                // Get position in AMReX convention to calculate corresponding index.
                // Ideally this will be replaced with the AMReX NGP interpolator
                // Always do x direction. No RZ case because it's not implemented, and code
                // will have aborted
                int ii = 0, jj = 0, kk = 0;
                amrex::ParticleReal x = p.pos(0);
                amrex::Real lx = (x - plo[0]) * dxi[0];
                ii = static_cast<int>(amrex::Math::floor(lx));
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D)
                amrex::ParticleReal y = p.pos(1);
                amrex::Real ly = (y - plo[1]) * dxi[1];
                jj = static_cast<int>(amrex::Math::floor(ly));
#endif
#if defined(WARPX_DIM_3D)
                amrex::ParticleReal z = p.pos(2);
                amrex::Real lz = (z - plo[2]) * dxi[2];
                kk = static_cast<int>(amrex::Math::floor(lz));
#endif

                // Fix dimensions since parser assumes u = gamma * v / c
                amrex::ParticleReal ux = p.rdata(PIdx::ux) / PhysConst::c;
                amrex::ParticleReal uy = p.rdata(PIdx::uy) / PhysConst::c;
                amrex::ParticleReal uz = p.rdata(PIdx::uz) / PhysConst::c;
                amrex::Real value = map_fn(xw, yw, zw, ux, uy, uz);
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 0), p.rdata(PIdx::w) * value);
            });
    // Add the weight for each particle -- total number of particles of this species
    ParticleToMesh(pc, ppc_mf, m_lev,
            [=] AMREX_GPU_DEVICE (const WarpXParticleContainer::SuperParticleType& p,
                amrex::Array4<amrex::Real> const& out_array,
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
            {

                // Get position in AMReX convention to calculate corresponding index.
                // Ideally this will be replaced with the AMReX NGP interpolator
                // Always do x direction. No RZ case because it's not implemented, and code
                // will have aborted
                int ii = 0, jj = 0, kk = 0;
                amrex::ParticleReal x = p.pos(0);
                amrex::Real lx = (x - plo[0]) * dxi[0];
                ii = static_cast<int>(amrex::Math::floor(lx));
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D)
                amrex::ParticleReal y = p.pos(1);
                amrex::Real ly = (y - plo[1]) * dxi[1];
                jj = static_cast<int>(amrex::Math::floor(ly));
#endif
#if defined(WARPX_DIM_3D)
                amrex::ParticleReal z = p.pos(2);
                amrex::Real lz = (z - plo[2]) * dxi[2];
                kk = static_cast<int>(amrex::Math::floor(lz));
#endif
                amrex::Gpu::Atomic::AddNoRet(&out_array(ii, jj, kk, 0), p.rdata(PIdx::w));
            });
    // Divide value by number of particles for average. Set average to zero if there are no particles
    for (amrex::MFIter mfi(red_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& box = mfi.tilebox();
        amrex::Array4<amrex::Real> const& a_red = red_mf.array(mfi);
        amrex::Array4<amrex::Real const> const& a_ppc = ppc_mf.const_array(mfi);
        amrex::ParallelFor(box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (a_ppc(i,j,k,0) == 0) a_red(i,j,k,0) = 0;
                    else a_red(i,j,k,0) = a_red(i,j,k,0) / a_ppc(i,j,k,0);
                });
    }

    // Coarsen and interpolate from ppc_mf to the output diagnostic MultiFab, mf_dst.
    CoarsenIO::Coarsen(mf_dst, red_mf, dcomp, 0, nComp(), 0, m_crse_ratio);
}
