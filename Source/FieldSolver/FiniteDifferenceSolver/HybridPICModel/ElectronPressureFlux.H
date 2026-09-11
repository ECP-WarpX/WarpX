/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Prabhat Kumar (Helion Energy)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_HYBRIDPIC_ELECTRONPRESSUREFLUX_H_
#define WARPX_HYBRIDPIC_ELECTRONPRESSUREFLUX_H_

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Dim3.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

/** Device-callable building blocks of the conservative-flux electron pressure advance
 *
 *      d_t pe = -gamma div(ue pe) + (gamma-1) [ue.grad(pe) + Q + div(kappa grad(pe/n))]
 *
 *  shared by the theta-implicit in-loop advance and the explicit sub-stepped fluid
 *  advance, so that both time integrators use identical spatial operators. */
struct HybridPeFluxParams
{
    /** enthalpy-flux face reconstruction: 0 central, 1 van Albada limited upwind */
    int adv = 0;
};

/** Keeps the van Albada slope average finite where both one-sided slopes vanish. */
constexpr amrex::Real HybridPeVanAlbadaEps = 1.0e-100;

/** pe at the face between the p0 and p1 nodes (pm1/p2 = next nodes outward).
 *  adv == 1: donor value plus the van Albada limited half-slope, reconstructed from the
 *  upwind side of the face velocity uf (C1 in the state). Otherwise the central face
 *  value 0.5 (p0 + p1). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real HybridPeFace (HybridPeFluxParams const& p, amrex::Real uf,
                          amrex::Real pm1, amrex::Real p0,
                          amrex::Real p1, amrex::Real p2)
{
    using namespace amrex::literals;
    if (p.adv == 1) {
        const bool wind_p = (uf >= 0._rt);
        const amrex::Real q1 = wind_p ? pm1 : p2;
        const amrex::Real q2 = wind_p ? p0  : p1;   // donor
        const amrex::Real q3 = wind_p ? p1  : p0;
        const amrex::Real dq1 = q2 - q1;
        const amrex::Real dq2 = q3 - q2;
        return q2 + 0.5_rt * dq1*dq2*(dq1 + dq2)
            / (dq1*dq1 + dq2*dq2 + HybridPeVanAlbadaEps);
    }
    return 0.5_rt*(p0 + p1);
}

/** Positive-definite dissipative deposit of the collocated pairing
 *      Q = eta(rho,|Jp|,t) |Jp|^2 + eta_H(rho,|B|) sum_c sum_d <(d_d Jp_c)^2>_faces
 *  with one-sided gradients at non-periodic domain faces. Jp, B, rho nodal. */
template <class EtaEx, class EtaHEx, class DxiArr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real HybridPeJouleQ (int i, int j, int k,
                            amrex::Array4<amrex::Real const> const& jpx,
                            amrex::Array4<amrex::Real const> const& jpy,
                            amrex::Array4<amrex::Real const> const& jpz,
                            amrex::Array4<amrex::Real const> const& bxa,
                            amrex::Array4<amrex::Real const> const& bya,
                            amrex::Array4<amrex::Real const> const& bza,
                            amrex::Array4<amrex::Real const> const& rr,
                            EtaEx const& eta_ex, EtaHEx const& etah_ex,
                            bool inc_hyp, amrex::Real t_now,
                            DxiArr const& dxiq,
                            amrex::Dim3 const& qlo, amrex::Dim3 const& qhi,
                            amrex::GpuArray<bool, 3> const& qper)
{
    using namespace amrex::literals;
    const amrex::Real jx0 = jpx(i,j,k);
    const amrex::Real jy0 = jpy(i,j,k);
    const amrex::Real jz0 = jpz(i,j,k);
    const amrex::Real j2 = jx0*jx0 + jy0*jy0 + jz0*jz0;
    amrex::Real Q = eta_ex(rr(i,j,k,0), std::sqrt(j2), t_now) * j2;
    if (inc_hyp) {
        const amrex::Real Bmag = std::sqrt(
            bxa(i,j,k)*bxa(i,j,k) + bya(i,j,k)*bya(i,j,k)
            + bza(i,j,k)*bza(i,j,k));
        amrex::Real g2 = 0.0_rt;
        // face-averaged squared gradients; one-sided at non-periodic
        // domain faces (stale ghosts there)
        auto acc = [&] (amrex::Array4<amrex::Real const> const& a) {
            const int ip = (i < qhi.x || qper[0]) ? i+1 : i;
            const int im = (i > qlo.x || qper[0]) ? i-1 : i;
            amrex::Real dp = (a(ip,j,k) - a(i,j,k)) * dxiq[0];
            amrex::Real dm = (a(i,j,k) - a(im,j,k)) * dxiq[0];
            g2 += 0.5_rt*(dp*dp + dm*dm);
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D)
            const int jp2 = (j < qhi.y || qper[1]) ? j+1 : j;
            const int jm  = (j > qlo.y || qper[1]) ? j-1 : j;
            dp = (a(i,jp2,k) - a(i,j,k)) * dxiq[1];
            dm = (a(i,j,k) - a(i,jm,k)) * dxiq[1];
            g2 += 0.5_rt*(dp*dp + dm*dm);
#endif
#if defined(WARPX_DIM_3D)
            const int kp = (k < qhi.z || qper[2]) ? k+1 : k;
            const int km = (k > qlo.z || qper[2]) ? k-1 : k;
            dp = (a(i,j,kp) - a(i,j,k)) * dxiq[2];
            dm = (a(i,j,k) - a(i,j,km)) * dxiq[2];
            g2 += 0.5_rt*(dp*dp + dm*dm);
#endif
        };
        acc(jpx); acc(jpy); acc(jpz);
        Q += etah_ex(rr(i,j,k,0), Bmag) * g2;
    }
    return Q;
}

/** C1 positivity floor with compact support: identity for pe >= 2 eps, Huber blend
 *  below, value eps and slope 0 at pe <= 0. Energy-neutral wherever pe >= 2 eps. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real HybridPeFloor (amrex::Real pe_new, amrex::Real pe_eps)
{
    using namespace amrex::literals;
    if (pe_new >= 2.0_rt*pe_eps) { return pe_new; }
    if (pe_new > 0.0_rt) { return pe_eps + 0.25_rt*pe_new*pe_new/pe_eps; }
    return pe_eps;
}

/** Electron pressure of the floored adiabat, n_floor Te0 (n_floor/n0_ref)^(gamma-1):
 *  the pressure the gamma-law closure assigns to a cell at the density floor
 *  (Te0 in J). Reference scale of the positivity floor and of the vacuum halo pin. */
[[nodiscard]] inline
amrex::Real HybridPeFlooredAdiabat (amrex::Real n_floor, amrex::Real Te0_J,
                                    amrex::Real n0_ref, amrex::Real gamma)
{
    using namespace amrex::literals;
    return n_floor * Te0_J * std::pow(n_floor / n0_ref, gamma - 1.0_rt);
}

#endif // WARPX_HYBRIDPIC_ELECTRONPRESSUREFLUX_H_
