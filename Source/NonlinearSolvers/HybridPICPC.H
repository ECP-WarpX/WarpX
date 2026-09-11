/* Copyright 2026 Prabhat Kumar
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_HYBRID_PIC_PC_H_
#define WARPX_HYBRID_PIC_PC_H_

#include "Preconditioner.H"

#include "EmbeddedBoundary/Enabled.H"
#include "Fields.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/utils/Enums.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

/**
 * \brief Physics-based EMHD preconditioner for the theta-implicit hybrid-PIC
 *  solver, after L. Chacon, JCP 530 (2025) 113910.
 *
 *  With particles, density and electron pressure frozen, the only dependence
 *  of E_Ohm on the Newton iterate E is through Faraday's law,
 *  dB = -theta*dt curl(dE), giving the stiff electron (EMHD) block
 *
 *    J_pc dE = dE + U (curl curl dE),
 *    U v     = c_w (v x bhat) + a_e v - a_H del2 v,
 *    c_w     = theta*dt |B0| / (mu0 rho_c),  a_e = theta*dt eta / mu0,
 *    a_H     = theta*dt eta_H / mu0,
 *  with eta frozen at the linearization state: eta(rho, |J0|),
 *  J0 = curl B0 / mu0 (freezing at |J| = 0 would underestimate anomalous,
 *  J-dependent resistivity precisely in current sheets).
 *
 *  Reduced W-system. Because the dE coefficient is the identity, dE is
 *  eliminated EXACTLY in favor of W = curl curl dE -- the parabolization of
 *  Chacon 2025 Sec. 3 applied globally:
 *
 *    (I + curlcurl U) W = curlcurl b,      dE = b - U W.
 *
 *  The curl-free (electrostatic) content of b -- physical and large in
 *  hybrid PIC -- passes through the reconstruction exactly instead of being
 *  1/beta-suppressed by the defect-corrected block diagonal (beta = l2 c_w,
 *  the whistler CFL at the grid Nyquist; l2 = sum_d 2/dx_d^2).
 *
 *  The W-system is smoothed with damped block-Jacobi (Chacon 2025 Sec. 4.1):
 *  the defect correction builds the per-node diagonal with curlcurl ~ -del2
 *  while all residuals use the exact central-difference operators, so the
 *  fixed point is exact. The diagonal P = alpha I + beta (. x bhat), alpha =
 *  1 + l2 (a_e + a_H l2), has the analytic inverse
 *
 *    P^-1 v = [alpha v - beta (v x bhat) + (beta^2/alpha)(v.bhat) bhat]
 *             / (alpha^2 + beta^2).
 *
 *  Primal pair (pc_hybrid_pic.pair, the default on the collocated grid).
 *  The frozen electron-convection term of the E-row Jacobian,
 *  C v = (theta*dt ue0) x curl(v), ue0 = (J_i - J_plasma)0/rho_c, is kept
 *  exactly by solving the coupled pair
 *
 *    E-row: xE + (I+S)^-1 [ U xW - theta*dt ue0 x (curl xE) ] = b~
 *    W-row: xW - curlcurl xE = 0
 *
 *  with a per-node joint 2x2 block damped-Jacobi smoother (defect curlcurl
 *  -> l2 in the block only; residuals exact):
 *
 *    P dxE = rE - U_loc rW / s_W,   dxW = rW + s_W l2 dxE,
 *    P     = alpha' I + beta (. x bhat),
 *    alpha' = 1 + theta*dt |ue0| l1 + l2 (a_e + 2 l2 a_H),  beta = wdf l2 c_w,
 *
 *  where l1 = sum_d 2/dx_d, the hyper-resistive defect uses the compact-
 *  Laplacian Nyquist symbol 2 l2, and wdf = pair_whistler_defect places the
 *  whistler coupling defect inside the wide central curl-curl symbol (which
 *  peaks at l2/2). The W-half is stored scaled by the local weight
 *  s_W = max(c_w, 1/l2) (coefficient comp 14), so the W-residual of each
 *  node is weighted by its own E-row coupling ||U_loc||: this makes the
 *  block norm E-space accurate on strong density-contrast profiles. The
 *  electrostatic content lives on the E-row identity with an O(1) diagonal,
 *  and the whistler stiffness is inverted analytically inside the block.
 *  The m_xF/bF/rF containers hold the E-half, m_xW/bW/rW the (scaled)
 *  W-half. Requires mode = nested on a Cartesian grid.
 *
 *  Modes (pc_hybrid_pic.mode):
 *
 *  "smoother" (default): nsweeps damped Jacobi sweeps on the fine grid of
 *    the reduced W-system. A fixed linear operator, valid under the plain
 *    (non-flexible) amrex_gmres outer solver.
 *
 *  "nested": geometric-MG V(3,3) (factor-2 coarsening while the boxes
 *    allow, same DistributionMapping; coefficient injection, l2 rebuilt per
 *    level; nodal full-weighting restriction, (tri)linear prolongation; 8
 *    bottom sweeps) preconditions an inner right-preconditioned GMRES
 *    (CGS2 orthogonalization) on the W-system or the pair. The Krylov space
 *    absorbs the per-level Nyquist rings of the central-difference curl that
 *    smoothing alone cannot. In the W-system the E-residual of the
 *    reconstructed solution is r_E = -U r_W, so the inner tolerance is
 *    divided by max(1, beta_max) to target E-space accuracy; the pair block
 *    norm is E-space accurate as is. The inner solve makes the
 *    preconditioner NONLINEAR: it requires a flexible outer Krylov method
 *    (newton.linear_solver = petsc_ksp with PETSC_OPTIONS="-ksp_type
 *    fgmres").
 *
 *  Physical boundaries: PEC(_Insulator) faces are modeled exactly (identity
 *  tangential-E rows with dead columns and parity-mirrored delta-B/delta-J
 *  ghosts, matching the residual's boundary application; see ClosureClass);
 *  other non-periodic BCs fall back to a homogeneous zero-extension closure.
 *  With embedded boundaries, EB-covered nodes are masked off the same nodal
 *  eb_update_E flag the residual's Ohm solves consume and act as identity
 *  rows.
 *
 *  With implicit_evolve.use_mass_matrices_pc = 1 the deposited reduced
 *  mass-matrix diagonal Sigma = dJ/dE (the ion response, as consumed by the
 *  EM-PIC preconditioners) is folded through the Hall term as
 *  J_pc = (I + S) + U K, S v = (Sigma v) x B0 / rho_c, with the exact
 *  reduction preserved via the pointwise closed-form (I+S)^-1.
 *
 *  Geometries: Cartesian 1D/2D/3D and RZ restricted to the m = 0 azimuthal
 *  mode (the residual's CylindricalYee up/down stencil compositions,
 *  on-axis branches and axis mirror parities are reproduced exactly; the
 *  r = 0 axis and PEC outer wall are part of the boundary closure), each on
 *  the collocated (nodal) or the staggered (Yee) grid. On the staggered
 *  grid the internal objects keep the single-MultiFab layout on a
 *  CELL-indexed container with per-component staggered interpretation (E/J
 *  components cell-centered along their own direction, B components nodal
 *  along theirs; every component's index range fits the cell box, and
 *  periodic FillBoundary produces the correct images for both
 *  interpretations). Index-space the Yee curls are the compact up/down
 *  differences (in RZ the cylindrical stencils verbatim), the Hall rotation
 *  carries the residual's node-interpolation chain (HybridPICPCUStag), and
 *  the compact curl-curl composition doubles the Nyquist symbol, so the
 *  defect diagonal uses 4/dx^2 per direction as in RZ. The staggered pair
 *  (periodic BCs only) uses the node-sandwich block of SmoothLevel.
 *
 *  Input parameters (prefix pc_hybrid_pic.):
 *    mode                  (string, "smoother")  "smoother" | "nested"
 *    pair                  (bool, nested && collocated && Cartesian)
 *                                                primal-pair formulation
 *    sigma                 (float, 1.3 pair / 0.55 otherwise)
 *                                                Jacobi damping, in (0, 2)
 *    sigma_w               (float, 1.2 pair / = sigma)
 *                                                damping of the pair W rows
 *    pair_whistler_defect  (float, 0.4)          whistler defect wdf, in (0, 1]
 *    mg_floor              (int, 4)              coarsest MG box width (nested)
 *    inner_max             (int, 8)              inner-GMRES iteration cap
 *    inner_rtol            (float, 1.0)          inner-GMRES relative tolerance
 *                                                (pair: in the s_W-weighted
 *                                                block norm; W-system: divided
 *                                                by max(1, beta_max)); 1.0
 *                                                stops the pair solve after its
 *                                                first iteration, i.e. one
 *                                                V-cycle plus one Krylov step
 *                                                per application
 *    nsweeps               (int, 48)             sweeps in smoother mode and
 *                                                on a non-coarsenable grid
 *    verbose               (bool, false)
 *  and implicit_evolve.pe_ue_cap_fac (float, 1.0): the residual's marker-CFL
 *  cap of |ue| in the pe advance, mirrored by the pe-row leg.
 */
/** Apply (I + S)^-1 to (vx, vy, vz) in place, where S u = (diag(w) u) x bhat
 *  is the deposited-mass-matrix ion response folded through the Hall term
 *  (w_c = Sigma_cc |B0| / rho_c; identity when w = 0). Closed-form 3x3
 *  adjugate inverse; reduces to the analytic rotation inverse
 *  (v - m v x bhat + m^2 (v.bhat) bhat)/(1 + m^2) when w is isotropic. */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void HybridPICPCInvIpS (amrex::Real bhx, amrex::Real bhy, amrex::Real bhz,
                        amrex::Real wx, amrex::Real wy, amrex::Real wz,
                        amrex::Real& vx, amrex::Real& vy, amrex::Real& vz)
{
    using namespace amrex::literals;
    // off-diagonal entries of A = I + S: (S u)_x = wy uy bhz - wz uz bhy, etc.
    const amrex::Real M01 =  wy*bhz;
    const amrex::Real M02 = -wz*bhy;
    const amrex::Real M10 = -wx*bhz;
    const amrex::Real M12 =  wz*bhx;
    const amrex::Real M20 =  wx*bhy;
    const amrex::Real M21 = -wy*bhx;
    const amrex::Real a00 = 1.0_rt  - M12*M21;
    const amrex::Real a01 = M02*M21 - M01;
    const amrex::Real a02 = M01*M12 - M02;
    const amrex::Real a10 = M12*M20 - M10;
    const amrex::Real a11 = 1.0_rt  - M02*M20;
    const amrex::Real a12 = M02*M10 - M12;
    const amrex::Real a20 = M10*M21 - M20;
    const amrex::Real a21 = M01*M20 - M21;
    const amrex::Real a22 = 1.0_rt  - M01*M10;
    const amrex::Real invdet = 1.0_rt / (a00 + M01*a10 + M02*a20);
    const amrex::Real ux = vx, uy = vy, uz = vz;
    vx = (a00*ux + a01*uy + a02*uz)*invdet;
    vy = (a10*ux + a11*uy + a12*uz)*invdet;
    vz = (a20*ux + a21*uy + a22*uz)*invdet;
}

/** pair mode (see HybridPICPC::m_pair): inline evaluation of
 *      tc = (I+S)^-1 ( U(aW)/sW - w x F )   and   kk = sW * curl(F),
 *  with F = curl xE staged in a 3-comp array, aW the scaled auxiliary
 *  W-hat, w = theta*dt ue0 (coef comps 11-13) and sW the local W-hat scale
 *  (coef comp 14). Collocated Cartesian central stencils; one call per node
 *  inside each consumer (smoother update, residual, matvec). */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void HybridPICPCPairTcKk (amrex::Array4<const amrex::Real> const& coef,
                          amrex::Array4<const amrex::Real> const& aW,
                          amrex::Array4<const amrex::Real> const& F,
                          int i, int j, int k,
                          bool ion,
                          amrex::GpuArray<amrex::Real,3> const& invdx2,
                          amrex::GpuArray<amrex::Real,3> const& inv2dx,
                          amrex::Real& tcx, amrex::Real& tcy, amrex::Real& tcz,
                          amrex::Real& kkx, amrex::Real& kky, amrex::Real& kkz)
{
    using namespace amrex::literals;

    const amrex::Real bhx = coef(i,j,k,0);
    const amrex::Real bhy = coef(i,j,k,1);
    const amrex::Real bhz = coef(i,j,k,2);
    const amrex::Real cw = coef(i,j,k,3);
    const amrex::Real ae = coef(i,j,k,4);
    const amrex::Real ah = coef(i,j,k,5);
    const amrex::Real sW = coef(i,j,k,14);

    const amrex::Real ax = aW(i,j,k,0);
    const amrex::Real ay = aW(i,j,k,1);
    const amrex::Real az = aW(i,j,k,2);

    amrex::Real lx = 0.0_rt, ly = 0.0_rt, lz = 0.0_rt;
#if defined(WARPX_DIM_3D)
    lx = (aW(i+1,j,k,0) - 2.0_rt*ax + aW(i-1,j,k,0))*invdx2[0]
       + (aW(i,j+1,k,0) - 2.0_rt*ax + aW(i,j-1,k,0))*invdx2[1]
       + (aW(i,j,k+1,0) - 2.0_rt*ax + aW(i,j,k-1,0))*invdx2[2];
    ly = (aW(i+1,j,k,1) - 2.0_rt*ay + aW(i-1,j,k,1))*invdx2[0]
       + (aW(i,j+1,k,1) - 2.0_rt*ay + aW(i,j-1,k,1))*invdx2[1]
       + (aW(i,j,k+1,1) - 2.0_rt*ay + aW(i,j,k-1,1))*invdx2[2];
    lz = (aW(i+1,j,k,2) - 2.0_rt*az + aW(i-1,j,k,2))*invdx2[0]
       + (aW(i,j+1,k,2) - 2.0_rt*az + aW(i,j-1,k,2))*invdx2[1]
       + (aW(i,j,k+1,2) - 2.0_rt*az + aW(i,j,k-1,2))*invdx2[2];
    kkx = sW*((F(i,j+1,k,2) - F(i,j-1,k,2))*inv2dx[1]
            - (F(i,j,k+1,1) - F(i,j,k-1,1))*inv2dx[2]);
    kky = sW*((F(i,j,k+1,0) - F(i,j,k-1,0))*inv2dx[2]
            - (F(i+1,j,k,2) - F(i-1,j,k,2))*inv2dx[0]);
    kkz = sW*((F(i+1,j,k,1) - F(i-1,j,k,1))*inv2dx[0]
            - (F(i,j+1,k,0) - F(i,j-1,k,0))*inv2dx[1]);
#elif defined(WARPX_DIM_XZ)
    lx = (aW(i+1,j,k,0) - 2.0_rt*ax + aW(i-1,j,k,0))*invdx2[0]
       + (aW(i,j+1,k,0) - 2.0_rt*ax + aW(i,j-1,k,0))*invdx2[2];
    ly = (aW(i+1,j,k,1) - 2.0_rt*ay + aW(i-1,j,k,1))*invdx2[0]
       + (aW(i,j+1,k,1) - 2.0_rt*ay + aW(i,j-1,k,1))*invdx2[2];
    lz = (aW(i+1,j,k,2) - 2.0_rt*az + aW(i-1,j,k,2))*invdx2[0]
       + (aW(i,j+1,k,2) - 2.0_rt*az + aW(i,j-1,k,2))*invdx2[2];
    kkx = sW*(-(F(i,j+1,k,1) - F(i,j-1,k,1))*inv2dx[2]);
    kky = sW*((F(i,j+1,k,0) - F(i,j-1,k,0))*inv2dx[2]
            - (F(i+1,j,k,2) - F(i-1,j,k,2))*inv2dx[0]);
    kkz = sW*((F(i+1,j,k,1) - F(i-1,j,k,1))*inv2dx[0]);
#else
    lx = (aW(i+1,j,k,0) - 2.0_rt*ax + aW(i-1,j,k,0))*invdx2[2];
    ly = (aW(i+1,j,k,1) - 2.0_rt*ay + aW(i-1,j,k,1))*invdx2[2];
    lz = (aW(i+1,j,k,2) - 2.0_rt*az + aW(i-1,j,k,2))*invdx2[2];
    kkx = sW*(-(F(i+1,j,k,1) - F(i-1,j,k,1))*inv2dx[2]);
    kky = sW*((F(i+1,j,k,0) - F(i-1,j,k,0))*inv2dx[2]);
    kkz = 0.0_rt;
#endif
    // U(aW)/sW - w x F, then the pointwise ion rotation
    tcx = (cw*(ay*bhz - az*bhy) + ae*ax - ah*lx)/sW;
    tcy = (cw*(az*bhx - ax*bhz) + ae*ay - ah*ly)/sW;
    tcz = (cw*(ax*bhy - ay*bhx) + ae*az - ah*lz)/sW;
    const amrex::Real ex = coef(i,j,k,11);
    const amrex::Real ey = coef(i,j,k,12);
    const amrex::Real ez = coef(i,j,k,13);
    tcx -= ey*F(i,j,k,2) - ez*F(i,j,k,1);
    tcy -= ez*F(i,j,k,0) - ex*F(i,j,k,2);
    tcz -= ex*F(i,j,k,1) - ey*F(i,j,k,0);
    if (ion) {
        HybridPICPCInvIpS(bhx, bhy, bhz,
                          coef(i,j,k,6), coef(i,j,k,7), coef(i,j,k,8),
                          tcx, tcy, tcz);
    }
}

/** Staggered (Yee) rotation legs of U at container index (i,j,k), mirroring
 *  the residual's Ohm structure: the Hall term is computed at NODES (each
 *  J-like component down-averaged along its cell-centered dimension), then
 *  up-averaged back to each E staggering; the resistive and hyper-resistive
 *  legs act at the component's own location with the node-stored
 *  coefficients up-averaged there. The caller passes the compact
 *  (index-space) Laplacian terms, already evaluated at the component
 *  locations, and applies the ion-response rotation afterwards. */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void HybridPICPCUStag (amrex::Array4<amrex::Real const> const& av,
                       amrex::Array4<amrex::Real const> const& coef,
                       int i, int j, int k,
                       amrex::Real lx, amrex::Real ly, amrex::Real lz,
                       amrex::Real& ux, amrex::Real& uy, amrex::Real& uz)
{
    using namespace amrex::literals;
    // c-component of c_w (v~ x bhat) at node (mi,mj,mk)
    auto ncross = [&] (int mi, int mj, int mk, int c) -> amrex::Real
    {
#if defined(WARPX_DIM_3D)
        const amrex::Real vx = 0.5_rt*(av(mi-1,mj,mk,0) + av(mi,mj,mk,0));
        const amrex::Real vy = 0.5_rt*(av(mi,mj-1,mk,1) + av(mi,mj,mk,1));
        const amrex::Real vz = 0.5_rt*(av(mi,mj,mk-1,2) + av(mi,mj,mk,2));
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        const amrex::Real vx = 0.5_rt*(av(mi-1,mj,mk,0) + av(mi,mj,mk,0));
        const amrex::Real vy = av(mi,mj,mk,1);
        const amrex::Real vz = 0.5_rt*(av(mi,mj-1,mk,2) + av(mi,mj,mk,2));
#else // WARPX_DIM_1D_Z
        const amrex::Real vx = av(mi,mj,mk,0);
        const amrex::Real vy = av(mi,mj,mk,1);
        const amrex::Real vz = 0.5_rt*(av(mi-1,mj,mk,2) + av(mi,mj,mk,2));
#endif
        // like the residual's enE stage, the node cross carries NO
        // density: c_w rho undoes the 1/rho folded into coef(3), and the
        // caller divides by rho interpolated to the component location
        const amrex::Real cwr = coef(mi,mj,mk,3)*coef(mi,mj,mk,9);
        if (c == 0) {
            return cwr*(vy*coef(mi,mj,mk,2) - vz*coef(mi,mj,mk,1));
        } else if (c == 1) {
            return cwr*(vz*coef(mi,mj,mk,0) - vx*coef(mi,mj,mk,2));
        }
        return cwr*(vx*coef(mi,mj,mk,1) - vy*coef(mi,mj,mk,0));
    };
    // guarded divide: zero density (physical ghosts under the
    // zero-extension closure) means a zero Hall term, not NaN
    auto odiv = [] (amrex::Real num, amrex::Real den) {
        return (den > 0.0_rt) ? num/den : 0.0_rt;
    };
#if defined(WARPX_DIM_3D)
    ux = odiv(ncross(i,j,k,0) + ncross(i+1,j,k,0),
              coef(i,j,k,9) + coef(i+1,j,k,9))
       + 0.5_rt*(coef(i,j,k,4) + coef(i+1,j,k,4))*av(i,j,k,0)
       - 0.5_rt*(coef(i,j,k,5) + coef(i+1,j,k,5))*lx;
    uy = odiv(ncross(i,j,k,1) + ncross(i,j+1,k,1),
              coef(i,j,k,9) + coef(i,j+1,k,9))
       + 0.5_rt*(coef(i,j,k,4) + coef(i,j+1,k,4))*av(i,j,k,1)
       - 0.5_rt*(coef(i,j,k,5) + coef(i,j+1,k,5))*ly;
    uz = odiv(ncross(i,j,k,2) + ncross(i,j,k+1,2),
              coef(i,j,k,9) + coef(i,j,k+1,9))
       + 0.5_rt*(coef(i,j,k,4) + coef(i,j,k+1,4))*av(i,j,k,2)
       - 0.5_rt*(coef(i,j,k,5) + coef(i,j,k+1,5))*lz;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    ux = odiv(ncross(i,j,k,0) + ncross(i+1,j,k,0),
              coef(i,j,k,9) + coef(i+1,j,k,9))
       + 0.5_rt*(coef(i,j,k,4) + coef(i+1,j,k,4))*av(i,j,k,0)
       - 0.5_rt*(coef(i,j,k,5) + coef(i+1,j,k,5))*lx;
    uy = odiv(ncross(i,j,k,1), coef(i,j,k,9))
       + coef(i,j,k,4)*av(i,j,k,1) - coef(i,j,k,5)*ly;
    uz = odiv(ncross(i,j,k,2) + ncross(i,j+1,k,2),
              coef(i,j,k,9) + coef(i,j+1,k,9))
       + 0.5_rt*(coef(i,j,k,4) + coef(i,j+1,k,4))*av(i,j,k,2)
       - 0.5_rt*(coef(i,j,k,5) + coef(i,j+1,k,5))*lz;
#else // WARPX_DIM_1D_Z
    ux = odiv(ncross(i,j,k,0), coef(i,j,k,9))
       + coef(i,j,k,4)*av(i,j,k,0) - coef(i,j,k,5)*lx;
    uy = odiv(ncross(i,j,k,1), coef(i,j,k,9))
       + coef(i,j,k,4)*av(i,j,k,1) - coef(i,j,k,5)*ly;
    uz = odiv(ncross(i,j,k,2) + ncross(i+1,j,k,2),
              coef(i,j,k,9) + coef(i+1,j,k,9))
       + 0.5_rt*(coef(i,j,k,4) + coef(i+1,j,k,4))*av(i,j,k,2)
       - 0.5_rt*(coef(i,j,k,5) + coef(i+1,j,k,5))*lz;
#endif
}

/** Staggered (Yee) variant of HybridPICPCPairTcKk: same contract, with the
 *  three staggered replacements. (a) U(aW) goes through the
 *  HybridPICPCUStag node chain (the compact index-space Laplacian is legal
 *  at each component's own location). (b) kk = sW curl(F) uses the compact
 *  downward (Ampere) differences, landing on the E/J staggering. (c) the
 *  convection term w x F mirrors the Hall cross: F (B-staggered) is
 *  down-averaged to nodes over its cell-centered dimensions, crossed with
 *  the node-built w = theta*dt ue0 (coef comps 11-13), and up-averaged
 *  back to each E staggering. The ion rotation keeps the node-collocation
 *  shortcut of the staggered W-path: node coefficients applied to the
 *  edge-centered triple (an upper-bound treatment, safe for a PC). */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void HybridPICPCPairTcKkStag (amrex::Array4<const amrex::Real> const& coef,
                              amrex::Array4<const amrex::Real> const& aW,
                              amrex::Array4<const amrex::Real> const& F,
                              int i, int j, int k,
                              bool ion,
                              amrex::GpuArray<amrex::Real,3> const& invdx2,
                              amrex::GpuArray<amrex::Real,3> const& inv2dx,
                              amrex::Real& tcx, amrex::Real& tcy, amrex::Real& tcz,
                              amrex::Real& kkx, amrex::Real& kky, amrex::Real& kkz)
{
    using namespace amrex::literals;

    const amrex::Real sW = coef(i,j,k,14);
    const amrex::Real ax = aW(i,j,k,0);
    const amrex::Real ay = aW(i,j,k,1);
    const amrex::Real az = aW(i,j,k,2);

    // index-space Laplacian at each component's own staggering (all
    // samples of a component are dx-separated in every direction)
    amrex::Real lx = 0.0_rt, ly = 0.0_rt, lz = 0.0_rt;
#if defined(WARPX_DIM_3D)
    lx = (aW(i+1,j,k,0) - 2.0_rt*ax + aW(i-1,j,k,0))*invdx2[0]
       + (aW(i,j+1,k,0) - 2.0_rt*ax + aW(i,j-1,k,0))*invdx2[1]
       + (aW(i,j,k+1,0) - 2.0_rt*ax + aW(i,j,k-1,0))*invdx2[2];
    ly = (aW(i+1,j,k,1) - 2.0_rt*ay + aW(i-1,j,k,1))*invdx2[0]
       + (aW(i,j+1,k,1) - 2.0_rt*ay + aW(i,j-1,k,1))*invdx2[1]
       + (aW(i,j,k+1,1) - 2.0_rt*ay + aW(i,j,k-1,1))*invdx2[2];
    lz = (aW(i+1,j,k,2) - 2.0_rt*az + aW(i-1,j,k,2))*invdx2[0]
       + (aW(i,j+1,k,2) - 2.0_rt*az + aW(i,j-1,k,2))*invdx2[1]
       + (aW(i,j,k+1,2) - 2.0_rt*az + aW(i,j,k-1,2))*invdx2[2];
    {
        const amrex::Real fx = 2.0_rt*inv2dx[0];
        const amrex::Real fy = 2.0_rt*inv2dx[1];
        const amrex::Real fz = 2.0_rt*inv2dx[2];
        kkx = sW*((F(i,j,k,2) - F(i,j-1,k,2))*fy
                - (F(i,j,k,1) - F(i,j,k-1,1))*fz);
        kky = sW*((F(i,j,k,0) - F(i,j,k-1,0))*fz
                - (F(i,j,k,2) - F(i-1,j,k,2))*fx);
        kkz = sW*((F(i,j,k,1) - F(i-1,j,k,1))*fx
                - (F(i,j,k,0) - F(i,j-1,k,0))*fy);
    }
#elif defined(WARPX_DIM_XZ)
    lx = (aW(i+1,j,k,0) - 2.0_rt*ax + aW(i-1,j,k,0))*invdx2[0]
       + (aW(i,j+1,k,0) - 2.0_rt*ax + aW(i,j-1,k,0))*invdx2[2];
    ly = (aW(i+1,j,k,1) - 2.0_rt*ay + aW(i-1,j,k,1))*invdx2[0]
       + (aW(i,j+1,k,1) - 2.0_rt*ay + aW(i,j-1,k,1))*invdx2[2];
    lz = (aW(i+1,j,k,2) - 2.0_rt*az + aW(i-1,j,k,2))*invdx2[0]
       + (aW(i,j+1,k,2) - 2.0_rt*az + aW(i,j-1,k,2))*invdx2[2];
    {
        const amrex::Real fx = 2.0_rt*inv2dx[0];
        const amrex::Real fz = 2.0_rt*inv2dx[2];
        kkx = sW*(-(F(i,j,k,1) - F(i,j-1,k,1))*fz);
        kky = sW*((F(i,j,k,0) - F(i,j-1,k,0))*fz
                - (F(i,j,k,2) - F(i-1,j,k,2))*fx);
        kkz = sW*((F(i,j,k,1) - F(i-1,j,k,1))*fx);
    }
#else
    lx = (aW(i+1,j,k,0) - 2.0_rt*ax + aW(i-1,j,k,0))*invdx2[2];
    ly = (aW(i+1,j,k,1) - 2.0_rt*ay + aW(i-1,j,k,1))*invdx2[2];
    lz = (aW(i+1,j,k,2) - 2.0_rt*az + aW(i-1,j,k,2))*invdx2[2];
    {
        const amrex::Real fz = 2.0_rt*inv2dx[2];
        kkx = sW*(-(F(i,j,k,1) - F(i-1,j,k,1))*fz);
        kky = sW*((F(i,j,k,0) - F(i-1,j,k,0))*fz);
        kkz = 0.0_rt;
    }
#endif

    // U(aW) through the staggered Hall/resistive node chain
    amrex::Real ux, uy, uz;
    HybridPICPCUStag(aW, coef, i, j, k, lx, ly, lz, ux, uy, uz);

    // c-component of w x F at node (mi,mj,mk): F down-averaged to the node
    // over its cell-centered dimensions (the bn* pattern of the coefficient
    // build), w read at the node
    auto ncrossF = [&] (int mi, int mj, int mk, int c) -> amrex::Real
    {
#if defined(WARPX_DIM_3D)
        const amrex::Real fx = 0.25_rt*(F(mi,mj-1,mk-1,0) + F(mi,mj,mk-1,0)
                                      + F(mi,mj-1,mk,0)   + F(mi,mj,mk,0));
        const amrex::Real fy = 0.25_rt*(F(mi-1,mj,mk-1,1) + F(mi,mj,mk-1,1)
                                      + F(mi-1,mj,mk,1)   + F(mi,mj,mk,1));
        const amrex::Real fz = 0.25_rt*(F(mi-1,mj-1,mk,2) + F(mi,mj-1,mk,2)
                                      + F(mi-1,mj,mk,2)   + F(mi,mj,mk,2));
#elif defined(WARPX_DIM_XZ)
        const amrex::Real fx = 0.5_rt*(F(mi,mj-1,mk,0) + F(mi,mj,mk,0));
        const amrex::Real fy = 0.25_rt*(F(mi-1,mj-1,mk,1) + F(mi,mj-1,mk,1)
                                      + F(mi-1,mj,mk,1)   + F(mi,mj,mk,1));
        const amrex::Real fz = 0.5_rt*(F(mi-1,mj,mk,2) + F(mi,mj,mk,2));
#else
        const amrex::Real fx = 0.5_rt*(F(mi-1,mj,mk,0) + F(mi,mj,mk,0));
        const amrex::Real fy = 0.5_rt*(F(mi-1,mj,mk,1) + F(mi,mj,mk,1));
        const amrex::Real fz = F(mi,mj,mk,2);
#endif
        const amrex::Real wx = coef(mi,mj,mk,11);
        const amrex::Real wy = coef(mi,mj,mk,12);
        const amrex::Real wz = coef(mi,mj,mk,13);
        if (c == 0) { return wy*fz - wz*fy; }
        if (c == 1) { return wz*fx - wx*fz; }
        return wx*fy - wy*fx;
    };
    // up-average the node cross back to each E staggering
#if defined(WARPX_DIM_3D)
    const amrex::Real cvx = 0.5_rt*(ncrossF(i,j,k,0) + ncrossF(i+1,j,k,0));
    const amrex::Real cvy = 0.5_rt*(ncrossF(i,j,k,1) + ncrossF(i,j+1,k,1));
    const amrex::Real cvz = 0.5_rt*(ncrossF(i,j,k,2) + ncrossF(i,j,k+1,2));
#elif defined(WARPX_DIM_XZ)
    const amrex::Real cvx = 0.5_rt*(ncrossF(i,j,k,0) + ncrossF(i+1,j,k,0));
    const amrex::Real cvy = ncrossF(i,j,k,1);
    const amrex::Real cvz = 0.5_rt*(ncrossF(i,j,k,2) + ncrossF(i,j+1,k,2));
#else
    const amrex::Real cvx = ncrossF(i,j,k,0);
    const amrex::Real cvy = ncrossF(i,j,k,1);
    const amrex::Real cvz = 0.5_rt*(ncrossF(i,j,k,2) + ncrossF(i+1,j,k,2));
#endif

    tcx = ux/sW - cvx;
    tcy = uy/sW - cvy;
    tcz = uz/sW - cvz;
    if (ion) {
        HybridPICPCInvIpS(coef(i,j,k,0), coef(i,j,k,1), coef(i,j,k,2),
                          coef(i,j,k,6), coef(i,j,k,7), coef(i,j,k,8),
                          tcx, tcy, tcz);
    }
}

#if defined(WARPX_DIM_RZ)
/** RZ (m = 0) stencil helpers mirroring the residual's discrete operators
 *  on the collocated grid: the Faraday curl uses the CylindricalYee
 *  "upward" differences and the Ampere curl the "downward" ones, each with
 *  the residual's on-axis branches, and r evaluated at the
 *  staggering-assumed positions the kernels use (nodal for r/theta rows,
 *  cell-centered for the z-row of the Faraday curl and the r-row of the
 *  hyper-resistive Laplacian). Components are (r, theta, z) = (0, 1, 2);
 *  index dims are (r, z) = (i, j). */
namespace hybridpcrz {

    /** True curl of an E-like field, Faraday (upward) stencils. */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void CurlF (amrex::Array4<amrex::Real const> const& a,
                amrex::Real& cr, amrex::Real& ct, amrex::Real& cz,
                int i, int j, amrex::Real rmin, amrex::Real dr,
                amrex::Real inv_dr, amrex::Real inv_dz)
    {
        using namespace amrex::literals;
        const amrex::Real rn = rmin + i*dr;             // nodal position
        const amrex::Real rc = rmin + (i + 0.5_rt)*dr;  // cell-centered
        // (curl a)_r = -d_z a_theta; the Faraday kernel zeroes B_r on axis
        cr = (rn == 0.0_rt) ? 0.0_rt
           : -inv_dz*(a(i,j+1,0,1) - a(i,j,0,1));
        // (curl a)_theta = d_z a_r - d_r a_z
        ct = inv_dz*(a(i,j+1,0,0) - a(i,j,0,0))
           - inv_dr*(a(i+1,j,0,2) - a(i,j,0,2));
        // (curl a)_z = (1/r) d_r (r a_theta), r at the cell center
        cz = (1.0_rt/rc)*inv_dr*((rc + 0.5_rt*dr)*a(i+1,j,0,1)
                               - (rc - 0.5_rt*dr)*a(i,j,0,1));
    }

    /** True curl of a B-like field, Ampere (downward) stencils with the
     *  residual's on-axis regularizations. */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void CurlA (amrex::Array4<amrex::Real const> const& a,
                amrex::Real& cr, amrex::Real& ct, amrex::Real& cz,
                int i, int j, amrex::Real rmin, amrex::Real dr,
                amrex::Real inv_dr, amrex::Real inv_dz)
    {
        using namespace amrex::literals;
        const amrex::Real rn = rmin + i*dr;
        // (curl a)_r = -d_z a_theta
        cr = -inv_dz*(a(i,j,0,1) - a(i,j-1,0,1));
        if (rn > 0.5_rt*dr) {
            // (curl a)_theta = d_z a_r - d_r a_z
            ct = inv_dz*(a(i,j,0,0) - a(i,j-1,0,0))
               - inv_dr*(a(i,j,0,2) - a(i-1,j,0,2));
            // (curl a)_z = (1/r) d_r (r a_theta), r at the node
            cz = (1.0_rt/rn)*inv_dr*((rn + 0.5_rt*dr)*a(i,j,0,1)
                                   - (rn - 0.5_rt*dr)*a(i-1,j,0,1));
        } else {
            ct = 0.0_rt;                        // J_theta = 0 on axis
            cz = 4.0_rt*a(i,j,0,1)*inv_dr;      // J_z axis regularization
        }
    }

    /** Component-wise hyper-resistive Laplacian of a J-like field, the
     *  residual's RZ forms (r-row at the cell-centered radius including the
     *  -1/r^2 curvature term; theta-row zero on axis; z-row switching to
     *  the plain second difference on axis). */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void LapJ (amrex::Array4<amrex::Real const> const& a,
               amrex::Real& lr, amrex::Real& lt, amrex::Real& lz,
               int i, int j, amrex::Real rmin, amrex::Real dr,
               amrex::Real inv_dr, amrex::Real inv_dz)
    {
        using namespace amrex::literals;
        const amrex::Real rn = rmin + i*dr;
        const amrex::Real rc = rmin + (i + 0.5_rt)*dr;
        const amrex::Real inv_dr2 = inv_dr*inv_dr;
        const amrex::Real inv_dz2 = inv_dz*inv_dz;
        auto dzz = [&] (int c) {
            return inv_dz2*(a(i,j-1,0,c) - 2.0_rt*a(i,j,0,c) + a(i,j+1,0,c));
        };
        auto dr_rdr_over_r = [&] (amrex::Real r, int c) {
            return (1.0_rt/r)*inv_dr2*((r + 0.5_rt*dr)*(a(i+1,j,0,c) - a(i,j,0,c))
                                     - (r - 0.5_rt*dr)*(a(i,j,0,c) - a(i-1,j,0,c)));
        };
        lr = dr_rdr_over_r(rc, 0) + dzz(0) - a(i,j,0,0)/(rc*rc);
        lt = (rn > 0.0_rt)
           ? dr_rdr_over_r(rn, 1) + dzz(1) - a(i,j,0,1)/(rn*rn)
           : 0.0_rt;
        lz = dzz(2);
        if (rn > 0.5_rt*dr) {
            lz += dr_rdr_over_r(rn, 2);
        } else {
            lz += inv_dr2*(a(i-1,j,0,2) - 2.0_rt*a(i,j,0,2) + a(i+1,j,0,2));
        }
    }

} // namespace hybridpcrz
#endif // WARPX_DIM_RZ

template <class T, class Ops>
class HybridPICPC : public Preconditioner<T,Ops>
{
    public:

        using RT = typename T::value_type;
        using MF = amrex::MultiFab;

        /** Largest number of Krylov coefficients passed by value to the
         *  block kernels (GpuArray kernel argument); bounds inner_max. */
        static constexpr int kMaxKry = 64;
        /** Thread-block size of the batched dot-product kernel. */
        static constexpr int kDotThreads = 256;
        /** V-cycle smoothing depths (Chacon 2025 Sec. 4: V(3,3), 8 bottom
         *  sweeps). The pre-smooth carries the x-is-zero initialization, so
         *  kNuPre >= 1 is required. */
        static constexpr int kNuPre = 3;
        static constexpr int kNuPost = 3;
        static constexpr int kNuBottom = 8;
        /** Iterate ghost layers: one exchange feeds U (1) + curl (1) + curl
         *  (1) per sweep. MG coarsening keeps boxes >= mg_floor >= 4 cells
         *  so periodic exchanges never wrap more than one period. */
        static constexpr int kNGrow = 3;
        /** Number of frozen coefficient components (see m_coef). */
        static constexpr int kNCoef = 15;

        HybridPICPC () = default;
        ~HybridPICPC () override = default;

        HybridPICPC(const HybridPICPC&) = delete;
        HybridPICPC& operator=(const HybridPICPC&) = delete;
        HybridPICPC(HybridPICPC&&) noexcept = delete;
        HybridPICPC& operator=(HybridPICPC&&) noexcept = delete;

        void Define (const T&, Ops*) override;

        /** Refreeze the coefficients (bhat, c_w, eta terms) from the current
         *  Newton linearization state and inject them down the hierarchy. */
        void Update (const T& a_U) override;

        /** Approximately solve J_pc x = b via the reduced W-system or the
         *  primal pair. */
        void Apply (T& a_x, const T& a_b) override;

        void printParameters () const override;

        [[nodiscard]] inline bool IsDefined () const override { return m_is_defined; }

        /* Implementation details; public because CUDA requires the enclosing
         * function of an extended device lambda to have public access. */

        /** tmp = (I+S)^-1 U a on growntilebox(2); tmp2 = curl tmp on
         *  growntilebox(1). a must be ghost-filled (kNGrow layers). After
         *  this, curl(tmp2) is available on the valid region without further
         *  exchanges. */
        void StageKU (const MF& a, int a_lev);

        /** pair staging (see m_pair): m_tmp = F = curl xE on growntilebox(2)
         *  from the ghost-filled E-half m_xF[a_lev]; tc and kk are then
         *  evaluated inline per node (HybridPICPCPairTcKk) by every consumer. */
        void StagePair (int a_lev);

        /** Damped block-Jacobi sweeps on level a_lev: per sweep, one ghost
         *  exchange of the iterate(s), the staging, then the update. */
        void SmoothLevel (int a_lev, int a_nsweeps, bool a_x_is_zero);

        /** W-system: rW = bW - xW - curl curl (U xW); pair: rE = bF - xE - tc,
         *  rW = bW - xW + kk, on the valid region of a_lev. */
        void ResidualLevel (int a_lev);

        /** Inner-GMRES matvec on the fine level, in block form (see m_kryB):
         *  a_blk (ghost-filled, nkb comps) -> a_out (nkb comps, valid).
         *  W-system: out = a + curl curl (U a); pair: out_W = W - kk,
         *  out_E = xE + tc. */
        void ApplyReducedOp (const MF& a_blk, MF& a_out);

        /** Block helpers of the inner GMRES (see m_kryB). SplitBlock reads
         *  the nkb comps of a_blk starting at a_scomp into m_bW[0] (comps
         *  0..2) and m_bF[0] (comps 3..5) on the valid region; MergeBlock
         *  writes m_xW[0]/m_xF[0] into comps 0..nkb-1 of a_blk. */
        void SplitBlock (const MF& a_blk, int a_scomp);
        void MergeBlock (MF& a_blk) const;

        /** Batched masked inner products of the block a_w with the first
         *  a_nv Arnoldi basis vectors of m_kryB: a_out[n] = <a_w, v_n> for
         *  n < a_nv and, when a_with_ww, a_out[a_nv] = <a_w, a_w>. One
         *  kernel per box (per-block partial sums, no atomics), ONE
         *  device->host copy, fixed-order host summation (deterministic),
         *  one vector MPI reduction. */
        void BatchDot (const MF& a_w, int a_nv, bool a_with_ww, RT* a_out);

        /** Fused block linear combination on the valid region:
         *  dst(comps a_dcomp..) = a_scale * ( [a_src] - sum_{n<a_nv} a_c[n] v_n ),
         *  with the a_src term present iff a_src != nullptr. One kernel. */
        void BlockAxpy (MF& a_dst, int a_dcomp, const MF* a_src, int a_nv,
                        const amrex::GpuArray<amrex::Real, kMaxKry>& a_c,
                        amrex::Real a_scale) const;

        /** Full-weighting restriction of rW[a_lev] (and rF for the pair)
         *  into bW[a_lev+1] (bF[a_lev+1]). */
        void RestrictLevel (int a_lev);

        /** (Tri)linear prolongation of xW[a_lev+1] (xF) added into xW[a_lev]
         *  (xF). */
        void ProlongLevel (int a_lev);

        /** Coefficient injection from level a_lev-1 to a_lev. */
        void RestrictCoef (int a_lev);

        /** Which discrete curl composition stage this is: the Faraday curl
         *  (E-like -> B-like, upward stencils in RZ) or the Ampere curl
         *  (B-like -> J-like, downward stencils with the on-axis branches
         *  in RZ). Cartesian uses the same central stencil for both. */
        enum class CurlKind { EtoB, BtoJ };

        /** out = curl a on the valid region grown by a_grow (a must have
         *  ghosts), reading a's components a_scomp..a_scomp+2. */
        void ComputeCurl (const MF& a, MF& out, int a_lev, int a_grow,
                          CurlKind a_kind, int a_scomp = 0) const;

        /** Final reconstruction helper: out = U a on the valid region
         *  grown by a_grow (staggered outputs need one ghost ring so the
         *  top nodal planes of the reconstruction can be read back). */
        void ApplyUValid (const MF& a, MF& out, int a_lev, int a_grow) const;

        /** Pointwise (I+S)^-1 on the valid region grown by a_grow. */
        void RotateValid (MF& a, int a_grow) const;
        void RotateBInPlace ();

        /** Index-space copy of one solver-vector component into component
         *  a_comp of a 3-component container (MultiFab::Copy requires
         *  matching ixTypes; the cell-indexed container's valid region is a
         *  subset of every staggered component's, so the reads are safe). */
        void CopyCompIn (const MF& a_src, MF& a_dst, int a_comp, bool a_add) const;

        /** Staggered reconstruction output: x_n = bE_n - t_n on each
         *  component's own staggered valid box (top nodal planes read the
         *  containers' ghost rings, which hold the periodic images /
         *  closure values). */
        void WriteXOut (const MF& a_bE, const MF& a_t, T& a_x) const;

        /** Physical-boundary closure classes for the vector objects in the
         *  PC pipeline, mirroring the residual's PEC handling exactly
         *  (see ApplyDomainClosure): at a PEC face with normal component n,
         *   EField : linearized E-like objects (the PC input b, the pair's
         *            E-half and the linearized Ohm outputs U(.)): tangential
         *            components zeroed ON the wall plane (ApplyEfieldBoundary
         *            is applied to the iterate AND to the Ohm output in the
         *            residual), ghosts mirror-filled tangential-odd /
         *            normal-even;
         *   BField : first-curl (delta-B-like) intermediates: normal
         *            component zeroed on the wall, ghosts tangential-even /
         *            normal-odd (ApplyBfieldBoundary after the Faraday
         *            stage);
         *   JField : delta-J-like objects (the W-system iterates W = K dE =
         *            -mu0 dJ_plasma / (theta dt)): ghosts tangential-odd /
         *            normal-even, NO wall write (CalculatePlasmaCurrent
         *            output receives no BC; its wall/ghost values follow
         *            from the BC-filled B);
         *   BCoef  : the coefficient-state total B (mirror ghosts with
         *            B parities, no wall write -- the copied field already
         *            satisfies its BCs);
         *   Coef   : the frozen coefficients (comps 0-2 = bhat with B-parity
         *            signs, the rest even; mirror ghosts, no wall write). */
        enum class ClosureClass { EField, BField, JField, BCoef, Coef };

        /** Apply the physical-boundary closure of class a_cls to a on level
         *  a_lev (PEC/PEC_Insulator faces only; other non-periodic BCs keep
         *  the zero-extension), then zero EB-covered nodes on the fine level
         *  when an EB mask is present. The closure indexes components
         *  RELATIVE to a_scomp (a_ncomp < 0: all remaining components). */
        void ApplyDomainClosure (MF& a, int a_lev, ClosureClass a_cls,
                                 int a_scomp = 0, int a_ncomp = -1) const;

        /** Ghost exchange of a preconditioner container on level a_lev.
         *  When the level consists of a single box spanning the domain the
         *  AMReX FillBoundary reduces to a periodic wrap of the ghost region
         *  onto the valid box, done here as ONE kernel with no host-device
         *  synchronization. Equals FillBoundary for periodically consistent
         *  data (the only case in which FillBoundary is itself well defined
         *  on a nodal periodic box): only ghost cells that overlap a
         *  periodic image of the valid box are written -- nothing at
         *  non-periodic faces (the closures handle those), the valid region
         *  including the duplicated nodal plane is untouched, and the
         *  source of a ghost is its wrap into the primary [lo, lo+N) range.
         *  Falls back to FillBoundary for multi-box layouts. */
        void FillGhosts (MF& a, int a_lev) const;

        /** Enforce the exact identity rows of the true Jacobian on a_out:
         *  wall-tangential components at PEC faces get a_out = a_in, and
         *  EB-covered nodes get a_out = a_in (the PC acts as the identity
         *  there; the covered residual rows are structurally zero). */
        void OverrideIdentityRows (T& a_out, const T& a_in) const;

    private:

        /** Recursive V(kNuPre, kNuPost) cycle for the correction on level
         *  a_lev; a non-coarsenable grid smooths nsweeps times instead. */
        void Vcycle (int a_lev, bool a_x_is_zero);

        /** One preconditioner application for the inner GMRES (one V-cycle);
         *  result in m_xW[0] (and m_xF[0] for the pair). */
        void ApplyMToW (const MF& a_blk, int a_scomp);

        /** Inner GMRES on the W-system or the pair (<= inner_max its or
         *  inner_rtol); result in m_xW[0] (m_xF[0]). Classical Gram-Schmidt
         *  with reorthogonalization (CGS2): two batched projection passes
         *  per iteration -- two host syncs -- instead of the j+2 sequential
         *  reductions of modified GS. */
        void InnerGMRES ();

        /** Inner-GMRES budget and relative tolerance (nested mode). */
        int m_inner_max = 8;
        RT m_inner_rtol = 1.0;
        /** Jacobi damping of the smoother (E rows in pair mode). */
        RT m_sigma = 0.55;
        /** Damping of the pair W-row update. The W eigenvalue of P^-1 A is 1
         *  on curl-null content, so sigma_w > 1 under-damps W junk while
         *  sigma > 1 speeds the whistler band. */
        RT m_sigma_w = 0.55;

        /** Ion-response leg active: the deposited reduced mass-matrix
         *  diagonal (raw dJ/dE) is available via
         *  implicit_evolve.use_mass_matrices_pc = 1 (set in Update()). */
        bool m_have_mm = false;

        /** pe-row leg: when the solver vector carries the pe scalar row
         *  (implicit_evolve.pe_newton_unknown), precondition it with a
         *  block-upper-triangular leg instead of the bare identity:
         *  x_p = D_p^-1 b_p with the pointwise diagonal D_p = 1 + coef(10),
         *  coef(10) = theta*dt*[gamma*sum_d |ue_d|/dx_d + conduction], then
         *  b_E <- b_E - grad(x_p/kappa)/rho_e before the E-system solve (the
         *  Ohm gradient coupling dF_E/dpe the identity leg left entirely to
         *  GMRES). |ue| is clamped like the residual's pe advance, ue_c =
         *  c tanh(|ue|/c), c = implicit_evolve.pe_ue_cap_fac * dx_min /
         *  (theta dt). Collocated Cartesian, in-loop pe closure only;
         *  inert for E-only vectors. */
        bool m_pe_row_active = false;
        std::unique_ptr<MF> m_pe_row_scr;   //!< ghosted dpe scratch (1 comp)
        amrex::Real m_pe_ue_cap_fac = 1.0;

        /** Primal-pair formulation (see the class docs). */
        bool m_pair = false;
        /** Whistler coupling defect wdf: d_perp = wdf*l2 in the block
         *  (beta = d_perp*c_w) and in the perp-b part of the W increment
         *  (par-b keeps l2). Collocated pair only. */
        RT m_pair_wdf = RT(0.4);
        /** Pair two-field state: per-level E unknown/RHS/residual. */
        std::vector<std::unique_ptr<MF>> m_xF;
        std::vector<std::unique_ptr<MF>> m_bF;
        std::vector<std::unique_ptr<MF>> m_rF;
        /** Ping-pong partners of m_xW/m_xF for the collocated pair smoother:
         *  its fused update reads the xW neighbours (hyper-resistive
         *  Laplacian inside HybridPICPCPairTcKk) of the iterate it writes,
         *  so an in-place update is a scheduling-dependent Jacobi/Gauss-
         *  Seidel mix. Each sweep writes the partner and the unique_ptrs
         *  are swapped (no copy), making the sweep a deterministic Jacobi
         *  iteration. */
        std::vector<std::unique_ptr<MF>> m_xW2;
        std::vector<std::unique_ptr<MF>> m_xF2;

        bool m_is_defined = false;
        bool m_nested = false;
        int m_mg_floor = 4;
        RT m_beta_max = 0.0;
        int m_nsweeps = 48;
        bool m_verbose = false;
        /** Per level: single-kernel ghost fill eligible (see FillGhosts). */
        std::vector<bool> m_fastfill_ok;
        /** True when the level solves need no pre-zeroed x: the collocated
         *  pair sweep writes x = 0 + sigma dE on the first sweep of a zero
         *  initial guess and never reads x before (ghosts are refilled from
         *  the written valid region), so the setVal kernels per level per
         *  V-cycle are skipped. */
        bool m_zero_free = false;

        Ops* m_ops = nullptr;

        int m_nlevels = 1;
        amrex::Vector<amrex::Geometry> m_geom;

        /** Per level: l2 = sum_d 2/dx_d^2 (the diagonal of -del2). */
        std::vector<RT> m_l2;
        /** Per level: l1 = sum_d 2/dx_d (Nyquist symbol of the first
         *  derivative; sizes the pair block's convection diagonal). */
        std::vector<RT> m_l1;
        std::vector<amrex::GpuArray<amrex::Real, 3>> m_inv2dx;
        std::vector<amrex::GpuArray<amrex::Real, 3>> m_invdx2;

        /** Per-node frozen coefficients (per level, kNCoef comps):
         *  comp 0-2  : bhat (unit total B),
         *  comp 3    : c_w = theta*dt |B| / (mu0 rho_e),
         *  comp 4    : a_e = theta*dt eta(rho, |J0|) / mu0,
         *  comp 5    : a_H = theta*dt eta_H / mu0,
         *  comp 6-8  : w_c = Sigma_cc |B| / rho_e (deposited-mass-matrix
         *              ion response through the Hall term; zero unless
         *              implicit_evolve.use_mass_matrices_pc = 1),
         *  comp 9    : rho_e (the floored charge density; used by the
         *              staggered Hall chain, which -- like the residual --
         *              crosses at the nodes WITHOUT the density and divides
         *              by rho interpolated to the component location),
         *  comp 10   : pe-row diagonal D_p - 1 (level 0 only),
         *  comp 11-13: theta*dt ue0 = theta*dt (J_i - J_plasma)/rho_e
         *              (pair convection),
         *  comp 14   : local W-hat scale s_W = max(c_w, 1/l2) (pair). */
        std::vector<std::unique_ptr<MF>> m_coef;

        /** Per level (all 3-component): W iterate (kNGrow ghosts), RHS (0),
         *  residual (1, for restriction), and the two stage scratches. */
        std::vector<std::unique_ptr<MF>> m_xW;
        std::vector<std::unique_ptr<MF>> m_bW;
        std::vector<std::unique_ptr<MF>> m_rW;
        std::vector<std::unique_ptr<MF>> m_tmp;
        std::vector<std::unique_ptr<MF>> m_tmp2;

        /** Fine level: ghosted copy of the PC input b (2 ghosts: one
         *  exchange feeds the curl-curl RHS assembly). */
        std::unique_ptr<MF> m_bE;

        /** Inner-GMRES storage (nested mode), in BLOCK form: every Krylov
         *  object carries m_nkb components -- the W-half in comps 0..2 and,
         *  for the pair, the E-half in comps 3..5 -- so each vector
         *  operation is ONE kernel, and the whole Arnoldi basis lives in a
         *  single MultiFab (basis vector j = comps [j*m_nkb, (j+1)*m_nkb)),
         *  so the Gram-Schmidt projections onto all j+1 basis vectors are a
         *  single batched kernel. */
        int m_nkb = 3;
        std::unique_ptr<MF> m_kryB;   //!< (inner_max+1)*nkb comps, 0 ghosts
        std::unique_ptr<MF> m_bKB;    //!< RHS block, 0 ghosts
        std::unique_ptr<MF> m_wB;     //!< A M v_j, 0 ghosts
        std::unique_ptr<MF> m_zB;     //!< M v_j, kNGrow ghosts (matvec input)
        /** Batched dot-product scratch: per-block partial sums, slot-major
         *  [slot][global block], one device->host copy per pass. */
        amrex::Gpu::DeviceVector<amrex::Real> m_dot_dev;
        amrex::Gpu::PinnedVector<amrex::Real> m_dot_host;
        std::vector<int> m_dot_nblocks;      //!< per local box
        std::vector<int> m_dot_box_offset;   //!< prefix sum of m_dot_nblocks
        int m_dot_nblocks_total = 0;

        /** Physical-BC closure state: per index dimension and side, whether
         *  the face is a PEC(_Insulator) wall the closure must model; the
         *  index-dim -> vector-component map; and the EB mask (the nodal
         *  eb_update_E flag on the collocated grid; nullptr without EB). */
        bool m_has_pec_wall = false;
        /** RZ: the r-lo face is the symmetry axis (rmin == 0). */
        bool m_has_axis = false;
        /** Staggered (Yee) grid: solver-vector components carry their
         *  native staggerings; internally everything is stored on the
         *  nodal-superset layout with per-component staggered
         *  interpretation. E_c is cell-centered along the index dimension
         *  that maps to component c (m_dim2comp) and nodal elsewhere. */
        bool m_staggered = false;
        bool m_wall_pec[AMREX_SPACEDIM][2] = {};
        /** Non-periodic faces under the zero-extension fallback: on the
         *  staggered grid the hi-side wall plane of components nodal along
         *  the face normal lives in the container's ghost region and is not
         *  modeled -- those rows are made exact identity rows instead. */
        bool m_wall_open[AMREX_SPACEDIM][2] = {};
        int m_dim2comp[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 1, 2)};
        const amrex::iMultiFab* m_eb_flag = nullptr;

};

template <class T, class Ops>
void HybridPICPC<T,Ops>::Define ( const T& a_U, Ops* a_ops )
{
    using namespace amrex;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined, "HybridPICPC::Define() called on defined object");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_U.IsDefined(), "HybridPICPC::Define() called with undefined solver vector");

    m_ops = a_ops;

    auto& warpx = WarpX::GetInstance();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        warpx.get_pointer_HybridPICModel() != nullptr,
        "pc_hybrid_pic requires the hybrid-PIC (Ohm's law) field solver");
#if defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    WARPX_ABORT_WITH_MESSAGE("pc_hybrid_pic is not implemented for 1D radial geometries");
#elif defined(WARPX_DIM_RZ)
    // RZ support is m = 0 only, with the axis at r = 0 (the residual's
    // CylindricalYee stencils and on-axis branches are mirrored exactly)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::n_rz_azimuthal_modes == 1,
        "pc_hybrid_pic in RZ supports the m = 0 azimuthal mode only");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(warpx.Geom(0).ProbLo(0) == 0.0,
        "pc_hybrid_pic in RZ requires the axis (r = 0) in the domain");
#endif
    m_staggered =
        (WarpX::grid_type != ablastr::utils::enums::GridType::Collocated);
    if (m_staggered) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            WarpX::grid_type == ablastr::utils::enums::GridType::Staggered,
            "pc_hybrid_pic supports the collocated and staggered (Yee) "
            "grids only");
    }

    // Physical-BC closure: detect PEC(_Insulator) faces so the PC's operator
    // reproduces the residual's boundary rows (identity wall-tangential-E
    // rows with dead columns; parity-mirrored ghosts for the delta-B /
    // delta-J stages). Other non-periodic BCs keep the zero-extension
    // closure (warned once).
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    m_dim2comp[0] = 0; m_dim2comp[1] = 2;
#elif defined(WARPX_DIM_1D_Z)
    m_dim2comp[0] = 2;
#endif
#if defined(WARPX_DIM_RZ)
    m_has_axis = true;
#endif
    {
        const auto& fbc_lo = a_ops->GetFieldBoundaryLo();
        const auto& fbc_hi = a_ops->GetFieldBoundaryHi();
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            if (warpx.Geom(0).isPeriodic(d)) { continue; }
            for (int side = 0; side < 2; side++) {
                const FieldBoundaryType bt =
                    (side == 0) ? fbc_lo[d] : fbc_hi[d];
                if (bt == FieldBoundaryType::PEC
                    || bt == FieldBoundaryType::PEC_Insulator) {
                    m_wall_pec[d][side] = true;
                    m_has_pec_wall = true;
                } else {
                    m_wall_open[d][side] = true;
                    ablastr::warn_manager::WMRecordWarning("HybridPICPC",
                        std::string("pc_hybrid_pic: non-periodic, non-PEC "
                        "field boundary detected; the preconditioner uses ")
                        + (m_staggered
                           ? "a homogeneous zero-extension ghost closure with "
                             "exact identity rows on the unmodeled hi-side "
                             "wall planes on that face (approximate boundary "
                             "rows)."
                           : "a homogeneous zero-extension closure on that "
                             "face (approximate boundary rows)."));
                }
            }
        }
    }

    const amrex::ParmParse pp("pc_hybrid_pic");
    std::string mode = "smoother";
    pp.query("mode", mode);
    if (mode == "nested") {
        m_nested = true;
        // The inner GMRES stops at an RHS-dependent iteration count, so the
        // preconditioner is nonlinear: classical (fixed-M) right-PC GMRES is
        // silently invalid. Enforce a flexible outer solver; the KSP type
        // itself is checked after KSPSetFromOptions in WarpX_PETSc.cpp.
        std::string nlsolver = "newton";
        const amrex::ParmParse pp_ie("implicit_evolve");
        pp_ie.query("nonlinear_solver", nlsolver);
        if (nlsolver == "newton") {
            std::string linsolver = "amrex_gmres";
            const amrex::ParmParse pp_newton("newton");
            pp_newton.query("linear_solver", linsolver);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(linsolver == "petsc_ksp",
                "pc_hybrid_pic.mode = nested is a nonlinear preconditioner and "
                "requires a flexible outer Krylov solver: set "
                "newton.linear_solver = petsc_ksp and "
                "PETSC_OPTIONS=\"-ksp_type fgmres\"");
        }
        amrex::Print() << "pc_hybrid_pic: nested mode (MG + inner GMRES) -- the "
                       << "preconditioner is nonlinear; the outer Krylov solver must be "
                       << "flexible (newton.linear_solver = petsc_ksp, "
                       << "PETSC_OPTIONS=\"-ksp_type fgmres\")\n";
    } else {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(mode == "smoother",
            "pc_hybrid_pic.mode must be 'smoother' or 'nested'");
    }
    pp.query("nsweeps", m_nsweeps);
    pp.query("verbose", m_verbose);
    pp.query("mg_floor", m_mg_floor);
    pp.query("inner_max", m_inner_max);
    {
        amrex::Real irtol = static_cast<amrex::Real>(m_inner_rtol);
        pp.query("inner_rtol", irtol);
        m_inner_rtol = static_cast<RT>(irtol);
    }
    // primal pair: the default wherever it is supported (nested mode on the
    // collocated Cartesian grid); the reduced W-system otherwise
#if defined(WARPX_DIM_RZ)
    m_pair = false;
#else
    m_pair = m_nested && !m_staggered;
#endif
    pp.query("pair", m_pair);
    m_sigma = m_pair ? RT(1.3) : RT(0.55);
    pp.query("sigma", m_sigma);
    m_sigma_w = m_pair ? RT(1.2) : m_sigma;
    {
        amrex::Real sgw = static_cast<amrex::Real>(m_sigma_w);
        pp.query("sigma_w", sgw);
        m_sigma_w = static_cast<RT>(sgw);
        amrex::Real wdf = static_cast<amrex::Real>(m_pair_wdf);
        pp.query("pair_whistler_defect", wdf);
        m_pair_wdf = static_cast<RT>(wdf);
    }
    {
        // mirror the residual's marker-CFL cap factor (ThetaImplicitHybrid)
        const amrex::ParmParse pp_impl("implicit_evolve");
        pp_impl.query("pe_ue_cap_fac", m_pe_ue_cap_fac);
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_sigma > 0.0 && m_sigma < 2.0,
        "pc_hybrid_pic.sigma must be in (0, 2)");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_sigma_w > RT(0.0) && m_sigma_w < RT(2.0),
        "pc_hybrid_pic.sigma_w must be in (0, 2)");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_pair_wdf > RT(0.0) && m_pair_wdf <= RT(1.0),
        "pc_hybrid_pic.pair_whistler_defect must be in (0, 1]");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_nsweeps > 0,
        "pc_hybrid_pic.nsweeps must be positive");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_inner_max > 0 && m_inner_rtol > RT(0.0),
        "pc_hybrid_pic.inner_max/inner_rtol must be positive");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_mg_floor >= 4,
        "pc_hybrid_pic.mg_floor must be >= 4");
#if defined(WARPX_DIM_RZ)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_pair,
        "pc_hybrid_pic.pair is Cartesian only (1D/2D/3D, collocated or "
        "staggered); cylindrical/spherical grids are not supported");
#endif
    if (m_pair) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_nested,
            "pc_hybrid_pic.pair requires mode = nested");
        // Physical boundaries (PEC walls, open zero-extension faces) ride
        // the same closure machinery as the reduced W-path: the E-half is
        // EField-classed at every fill (dead tangential columns at PEC
        // faces), the RHS b~ is closured before its copy into the E-half of
        // m_bKB, and the output passes through OverrideIdentityRows. The
        // staggered pair's two-field staging and output exchange assume
        // periodic ghosts.
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !m_staggered || warpx.Geom(0).isAllPeriodic(),
            "pc_hybrid_pic.pair on the staggered grid requires fully "
            "periodic BCs (physical boundaries: collocated only)");
    }

    // Storage on the solver-vector layout (single AMR level); in nested mode
    // the MG hierarchy coarsens that grid while every box stays
    // coarsenable, keeping the DistributionMapping so transfers are local.
    const amrex::MultiFab* E0 = a_U.getArrayVec()[0][0];
    const amrex::DistributionMapping& dm = E0->DistributionMap();

    // Container layout: nodal on the collocated grid (the solver vector's
    // own staggering); CELL-indexed on the staggered grid (see class docs)
    amrex::Vector<amrex::BoxArray> ba_node;
    ba_node.push_back(amrex::convert(E0->boxArray(),
        m_staggered ? amrex::IntVect::TheCellVector()
                    : amrex::IntVect::TheNodeVector()));
    m_geom.clear();
    m_geom.push_back(warpx.Geom(0));

    if (m_nested) {
        amrex::BoxArray ba_cell = amrex::convert(ba_node[0], amrex::IntVect::TheCellVector());
        while (true) {
            const amrex::Geometry& fg = m_geom.back();
            // min box width mg_floor >= 4: periodic exchanges of kNGrow
            // ghosts must not wrap more than one period
            if (!fg.Domain().coarsenable(2, m_mg_floor)
                || !ba_cell.coarsenable(2, m_mg_floor)) { break; }
            ba_cell.coarsen(2);
            ba_node.push_back(amrex::convert(ba_cell, ba_node[0].ixType()));
            amrex::Array<int,AMREX_SPACEDIM> is_per{};
            for (int d = 0; d < AMREX_SPACEDIM; d++) { is_per[d] = fg.isPeriodic(d); }
            m_geom.emplace_back(amrex::coarsen(fg.Domain(), 2), fg.ProbDomain(),
                                static_cast<int>(fg.Coord()), is_per);
        }
        if (m_geom.size() == 1) {
            amrex::Print() << "pc_hybrid_pic: grid is not coarsenable; the inner GMRES "
                           << "will be preconditioned by " << m_nsweeps << " sweeps\n";
        }
    }
    m_nlevels = static_cast<int>(m_geom.size());

    // single-kernel ghost fill eligibility (FillGhosts): one box spanning
    // the whole (index-type-converted) domain on that level
    m_fastfill_ok.assign(m_nlevels, false);
    for (int lev = 0; lev < m_nlevels; lev++) {
        const amrex::BoxArray& ba = ba_node[lev];
        const amrex::Box dom_t = amrex::convert(m_geom[lev].Domain(), ba.ixType());
        m_fastfill_ok[lev] = (ba.size() == 1) && (ba[0] == dom_t);
    }
#if defined(WARPX_DIM_RZ)
    m_zero_free = false;
#else
    m_zero_free = m_pair && !m_staggered;   // EB (set in Update) checked at use
#endif

#if defined(WARPX_DIM_3D)
    const int dirs[3] = {0, 1, 2};
    constexpr int ndirs = 3;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const int dirs[2] = {0, 2};
    constexpr int ndirs = 2;
#else // WARPX_DIM_1D_Z
    const int dirs[1] = {2};
    constexpr int ndirs = 1;
#endif
    m_l2.assign(m_nlevels, 0.0);
    m_l1.assign(m_nlevels, 0.0);
    m_inv2dx.assign(m_nlevels, {0.0, 0.0, 0.0});
    m_invdx2.assign(m_nlevels, {0.0, 0.0, 0.0});
    for (int lev = 0; lev < m_nlevels; lev++) {
        const auto* dx = m_geom[lev].CellSize();
        for (int d = 0; d < ndirs; d++) {
            const amrex::Real h = dx[d];
            m_l1[lev] += 2.0_rt / h;
            m_inv2dx[lev][dirs[d]] = 0.5_rt / h;
            m_invdx2[lev][dirs[d]] = 1.0_rt / (h*h);
#if defined(WARPX_DIM_RZ)
            // the RZ curl-curl is the compact down/up composition, whose
            // Nyquist symbol reaches 4/h^2 per direction (the Cartesian
            // wide central composition peaks at half that): size the
            // defect diagonal accordingly so the sigma < 2 damping bound
            // keeps its Cartesian semantics. The same applies to the
            // staggered (Yee) composition in Cartesian geometry.
            m_l2[lev] += 4.0_rt / (h*h);
#else
            m_l2[lev] += (m_staggered ? 4.0_rt : 2.0_rt) / (h*h);
#endif
        }
    }

    m_coef.resize(m_nlevels);
    m_xW.resize(m_nlevels);
    m_bW.resize(m_nlevels);
    m_rW.resize(m_nlevels);
    m_tmp.resize(m_nlevels);
    m_tmp2.resize(m_nlevels);
    m_xF.resize(m_nlevels);
    m_bF.resize(m_nlevels);
    m_rF.resize(m_nlevels);
    if (m_pair && !m_staggered) {
        m_xW2.resize(m_nlevels);
        m_xF2.resize(m_nlevels);
    }

    for (int lev = 0; lev < m_nlevels; lev++) {
        const amrex::BoxArray& ba = ba_node[lev];
        // kNGrow coefficient ghost layers: the U-stage evaluates U in a
        // 2-ring around the valid region, and the staggered up-averages
        // read the coefficients one ring beyond that; internal ghosts are
        // exchanged in Update(), physical ghosts follow the closures
        m_coef[lev] = std::make_unique<MF>(ba, dm, kNCoef, kNGrow);
        m_coef[lev]->setVal(0.0);
        m_xW[lev]   = std::make_unique<MF>(ba, dm, 3, kNGrow);
        m_rW[lev]   = std::make_unique<MF>(ba, dm, 3, 1);
        m_tmp[lev]  = std::make_unique<MF>(ba, dm, 3, kNGrow);
        m_tmp2[lev] = std::make_unique<MF>(ba, dm, 3, kNGrow);
        m_bW[lev]   = std::make_unique<MF>(ba, dm, 3, 0);
        // Domain (non-periodic) ghosts stay zero: homogeneous Dirichlet
        // closure on every level (exchanges fill periodic/internal only)
        m_xW[lev]->setVal(0.0);
        m_rW[lev]->setVal(0.0);
        m_tmp[lev]->setVal(0.0);
        m_tmp2[lev]->setVal(0.0);
        m_bW[lev]->setVal(0.0);
        if (m_pair) {
            m_xF[lev] = std::make_unique<MF>(ba, dm, 3, kNGrow);
            m_bF[lev] = std::make_unique<MF>(ba, dm, 3, 0);
            m_rF[lev] = std::make_unique<MF>(ba, dm, 3, 1);
            m_xF[lev]->setVal(0.0);
            m_bF[lev]->setVal(0.0);
            m_rF[lev]->setVal(0.0);
        }
        if (m_pair && !m_staggered) {
            // ping-pong partners of the pair iterates (see m_xW2)
            m_xW2[lev] = std::make_unique<MF>(ba, dm, 3, kNGrow);
            m_xF2[lev] = std::make_unique<MF>(ba, dm, 3, kNGrow);
            m_xW2[lev]->setVal(0.0);
            m_xF2[lev]->setVal(0.0);
        }
    }
    m_bE = std::make_unique<MF>(ba_node[0], dm, 3, 2);
    m_bE->setVal(0.0);
    m_pe_row_scr = std::make_unique<MF>(ba_node[0], dm, 1, 1);
    m_pe_row_scr->setVal(0.0);

    if (m_nested) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_inner_max + 1 <= kMaxKry,
            "pc_hybrid_pic.inner_max must be below 64 (Krylov coefficients "
            "are passed to the block kernels by value)");
        // block Krylov storage (see m_kryB): W-half comps 0..2, E-half
        // comps 3..5 for the pair
        m_nkb = m_pair ? 6 : 3;
        m_kryB = std::make_unique<MF>(ba_node[0], dm, (m_inner_max + 1) * m_nkb, 0);
        m_bKB  = std::make_unique<MF>(ba_node[0], dm, m_nkb, 0);
        m_wB   = std::make_unique<MF>(ba_node[0], dm, m_nkb, 0);
        m_zB   = std::make_unique<MF>(ba_node[0], dm, m_nkb, kNGrow);
        m_kryB->setVal(0.0);
        m_bKB->setVal(0.0);
        m_wB->setVal(0.0);
        m_zB->setVal(0.0);
        // batched-dot launch layout: ~2 cells per thread, at most 512
        // blocks per box; partial sums are slot-major over all local boxes
        m_dot_nblocks.clear();
        m_dot_box_offset.clear();
        m_dot_nblocks_total = 0;
        for (amrex::MFIter mfi(*m_wB, false); mfi.isValid(); ++mfi) {
            const amrex::Long npts = mfi.validbox().numPts();
            int nb = static_cast<int>((npts + 2*kDotThreads - 1) / (2*kDotThreads));
            nb = std::max(1, std::min(nb, 512));
            m_dot_box_offset.push_back(m_dot_nblocks_total);
            m_dot_nblocks.push_back(nb);
            m_dot_nblocks_total += nb;
        }
        const std::size_t ndot = static_cast<std::size_t>(m_inner_max + 2)
            * static_cast<std::size_t>(std::max(m_dot_nblocks_total, 1));
        m_dot_dev.resize(ndot);
        m_dot_host.resize(ndot);
    }

    m_is_defined = true;
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::printParameters () const
{
    using namespace amrex;
    Print() << "Hybrid-PIC EMHD preconditioner (pc_hybrid_pic):\n";
    Print() << "  grid: " << (m_staggered ? "staggered (Yee)" : "collocated")
            << "\n";
    if (m_nested) {
        Print() << "  mode: nested -- MG V(" << kNuPre << "," << kNuPost << ") ("
                << kNuBottom << " bottom sweeps, coarsest "
                << m_geom.back().Domain().length(0) << ") inside GMRES(<= "
                << m_inner_max << " its, rtol " << m_inner_rtol;
        if (m_pair) {
            Print() << " in the s_W-weighted pair block norm";
        } else {
            Print() << "/max(1, beta_max) in W-space";
        }
        Print() << "); requires a flexible outer solver\n";
        Print() << "  levels: " << m_nlevels << "\n";
        for (int lev = 0; lev < m_nlevels; lev++) {
            Print() << "    level " << lev << ": " << m_geom[lev].Domain()
                    << ", diag l2 = " << m_l2[lev] << "\n";
        }
    } else {
        Print() << "  mode: smoother, nsweeps = " << m_nsweeps
                << " (scale ~1.5-3x the peak whistler CFL beta = l2 c_w)\n";
        Print() << "  diag l2: " << m_l2[0] << "\n";
    }
    if (m_pair) {
        Print() << "  formulation: primal pair (xE, W), joint block Jacobi "
                << "(convection exact)\n";
        Print() << "  damping: " << m_sigma << " (E rows), " << m_sigma_w
                << " (W rows, sigma_w)\n";
        Print() << "  pair block defects: hyper-resistive Nyquist 2*l2, whistler "
                << m_pair_wdf << "*l2 (pair_whistler_defect)\n";
    } else {
        Print() << "  formulation: reduced W-system\n";
        Print() << "  damping: " << m_sigma << "\n";
    }
    Print() << "  ion response (deposited mass-matrix diagonal, "
            << "implicit_evolve.use_mass_matrices_pc): "
            << ((m_ops != nullptr && m_ops->GetMassMatricesCoeff() != nullptr)
                ? "true" : "false")
            << "\n";
    Print() << "  boundary closure: "
            << (m_has_pec_wall
                ? "PEC wall model (identity tangential-E rows, dead "
                  "columns, parity-mirrored delta-B/delta-J ghosts)"
                : "periodic/zero-extension")
            << (EB::enabled() ? " + EB staircase mask (eb_update_E)" : "")
            << "\n";
    Print() << "  single-kernel ghost fill: "
            << (m_fastfill_ok.empty() || !m_fastfill_ok[0] ? "off" : "on")
            << " (single domain-spanning box per level)\n";
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::Update ( const T& a_U )
{
    using namespace amrex;
    using warpx::fields::FieldType;
    using ablastr::fields::Direction;

    BL_PROFILE("HybridPICPC::Update()");

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(), "HybridPICPC::Update() called on undefined object");
    amrex::ignore_unused(a_U);

    auto& warpx = WarpX::GetInstance();
    const HybridPICModel* hybrid = warpx.get_pointer_HybridPICModel();

    // EB mask (staircase): the nodal eb_update_E flag the residual's Ohm
    // solves consume; on the collocated grid all three components share the
    // nodal staggering, so component 0 serves as the per-node mask
    m_eb_flag = nullptr;
    if (EB::enabled()) {
        m_eb_flag = warpx.GetEBUpdateEFlag()[0][0].get();
    }

    // Bfield_fp holds the internal B^{n+theta} at this point of the Newton
    // iteration (external fields are subtracted at the end of ComputeRHS);
    // the whistler rides on the total field, so add B_external here.
    const bool add_external = hybrid->m_add_external_fields;
    ablastr::fields::VectorField B_fp = warpx.m_fields.get_alldirs(FieldType::Bfield_fp, 0);
    ablastr::fields::VectorField B_ext;
    if (add_external) {
        B_ext = warpx.m_fields.get_alldirs(FieldType::hybrid_B_fp_external, 0);
    }
    const amrex::MultiFab* rho = warpx.m_fields.get(FieldType::rho_fp, 0);

    // pe-row leg (see m_pe_row_active): Update() cannot see whether the
    // solver vector carries the scalar row, so activate on the model gates
    // only; Apply() additionally gates on the scalar's presence, making
    // this inert for E-only vectors.
#if defined(WARPX_DIM_RZ)
    m_pe_row_active = false;
#else
    m_pe_row_active = !m_staggered && !hybrid->m_implicit_use_algebraic_closure;
#endif
    const bool pe_row_on = m_pe_row_active;
    // the pair's convection term needs the frozen electron velocity
    const bool ue_on = m_pair;
    ablastr::fields::VectorField Ji_fp, Jp_fp;
    if (pe_row_on || ue_on) {
        Ji_fp = warpx.m_fields.get_alldirs(FieldType::current_fp, 0);
        Jp_fp = warpx.m_fields.get_alldirs(FieldType::hybrid_current_fp_plasma, 0);
    }
    const RT gamma_e = static_cast<RT>(hybrid->m_gamma);

    const RT theta_dt = m_ops->GetTheta() * this->m_dt;
    const RT t_now = this->m_time;
    const RT rho_floor = hybrid->m_n_floor * PhysConst::q_e;
    // same C1 smooth floor as the Ohm kernels: the PC coefficients see the
    // identical rho_lim field the residual differentiates
    const RT floor_w = static_cast<RT>(hybrid->m_n_floor_smooth_width) * rho_floor;
    const RT one_over_mu0 = 1.0_rt / PhysConst::mu0;
    const auto eta = hybrid->m_eta;
    const auto eta_h = hybrid->m_eta_h;
    const bool include_hyper = hybrid->m_include_hyper_resistivity_term;

    // J0 = curl(B_tot)/mu0 at the linearization state, for freezing the
    // (possibly J-dependent) resistivity at |J0| rather than |J| = 0
    // (m_bE doubles as ghosted scratch here; Apply() rewrites it)
    for (int n = 0; n < 3; n++) {
        if (m_staggered) {
            CopyCompIn(*B_fp[n], *m_bE, n, false);
            if (add_external) { CopyCompIn(*B_ext[n], *m_bE, n, true); }
        } else {
            amrex::MultiFab::Copy(*m_bE, *B_fp[n], 0, n, 1, 0);
            if (add_external) {
                amrex::MultiFab::Add(*m_bE, *B_ext[n], 0, n, 1, 0);
            }
        }
    }
    FillGhosts(*m_bE, 0);
    ApplyDomainClosure(*m_bE, 0, ClosureClass::BCoef);
    ComputeCurl(*m_bE, *m_tmp[0], 0, 0, CurlKind::BtoJ);   // m_tmp[0] = mu0 J0

    // Reduced mass-matrix diagonal (raw dJ/dE, deposited/synced each Newton
    // iteration in PreLinearSolve) for the optional ion-response leg,
    // provided through implicit_evolve.use_mass_matrices_pc = 1
    const auto* mmpc = m_ops->GetMassMatricesCoeff();
    m_have_mm = (mmpc != nullptr) && !mmpc->empty();
    const bool have_mm = m_have_mm;
    const amrex::MultiFab* Sxx_mf = have_mm ? (*mmpc)[0][0] : nullptr;
    const amrex::MultiFab* Syy_mf = have_mm ? (*mmpc)[0][1] : nullptr;
    const amrex::MultiFab* Szz_mf = have_mm ? (*mmpc)[0][2] : nullptr;
    const int cxx = have_mm ? (Sxx_mf->nComp() - 1)/2 : 0;
    const int cyy = have_mm ? (Syy_mf->nComp() - 1)/2 : 0;
    const int czz = have_mm ? (Szz_mf->nComp() - 1)/2 : 0;
    if (have_mm) {
        // the deposited diagonals are read index-space (per-component E
        // staggerings on the Yee grid, nodal when collocated): the grids
        // must align cell-wise with the coefficient container
        AMREX_ALWAYS_ASSERT(
            amrex::convert(Sxx_mf->boxArray(), amrex::IntVect::TheCellVector())
                == amrex::convert(m_coef[0]->boxArray(),
                                  amrex::IntVect::TheCellVector())
            && Sxx_mf->DistributionMap() == m_coef[0]->DistributionMap());
    }
    const bool stag = m_staggered;
    const auto dxi_pc = m_geom[0].InvCellSizeArray();
#if defined(WARPX_DIM_RZ)
    amrex::ignore_unused(dxi_pc);   // the pe-row leg is Cartesian only
#endif
    const RT l2_pc = m_l2[0];
    // marker-CFL clamp of |ue| for the pe-row diagonal (see m_pe_row_active):
    // c = pe_ue_cap_fac * dx_min / (theta dt), the residual's cap; <= 0 = none
    amrex::Real dx_min_pc = m_geom[0].CellSize(0);
    for (int d = 1; d < AMREX_SPACEDIM; ++d) {
        dx_min_pc = std::min(dx_min_pc, m_geom[0].CellSize(d));
    }
    const amrex::Real ue_cap_pc =
        static_cast<amrex::Real>(m_pe_ue_cap_fac * dx_min_pc / theta_dt);
    // representative constant conductivity for the PC leg: a spatially
    // varying kappa_e(rho,Te) is sampled at the ambient reference state
    // (n0_ref, Te0) -- the PC need not track the exact coefficient field
    const RT kappa_pc = static_cast<RT>(
        hybrid->m_has_kappa_e_expression
            ? hybrid->m_kappa(hybrid->m_n0_ref * PhysConst::q_e,
                              hybrid->m_elec_temp / PhysConst::q_e)
            : hybrid->m_kappa_e);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_coef[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& coef = m_coef[0]->array(mfi);
        Array4<const Real> const& Bx = B_fp[0]->const_array(mfi);
        Array4<const Real> const& By = B_fp[1]->const_array(mfi);
        Array4<const Real> const& Bz = B_fp[2]->const_array(mfi);
        Array4<const Real> const& r = rho->const_array(mfi);
        Array4<const Real> Bxe, Bye, Bze;
        if (add_external) {
            Bxe = B_ext[0]->const_array(mfi);
            Bye = B_ext[1]->const_array(mfi);
            Bze = B_ext[2]->const_array(mfi);
        }
        Array4<const Real> const& cJ = m_tmp[0]->const_array(mfi);
        Array4<const Real> Sxx, Syy, Szz;
        if (have_mm) {
            Sxx = Sxx_mf->const_array(mfi);
            Syy = Syy_mf->const_array(mfi);
            Szz = Szz_mf->const_array(mfi);
        }
        Array4<const Real> Jix, Jiy, Jiz, Jpx, Jpy, Jpz;
        if (pe_row_on || ue_on) {
            Jix = Ji_fp[0]->const_array(mfi);
            Jiy = Ji_fp[1]->const_array(mfi);
            Jiz = Ji_fp[2]->const_array(mfi);
            Jpx = Jp_fp[0]->const_array(mfi);
            Jpy = Jp_fp[1]->const_array(mfi);
            Jpz = Jp_fp[2]->const_array(mfi);
        }
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // staggered: down-average each B component to the node over
                // its cell-centered dimensions (the residual's ghost-filled
                // B provides the lo-side reads; on the axis the odd
                // B_theta ghosts average to the physical zero)
                auto bnx = [&] (Array4<const Real> const& B) {
#if defined(WARPX_DIM_3D)
                    return 0.25_rt*(B(i,j-1,k-1) + B(i,j,k-1)
                                  + B(i,j-1,k) + B(i,j,k));
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    return 0.5_rt*(B(i,j-1,k) + B(i,j,k));
#else
                    return 0.5_rt*(B(i-1,j,k) + B(i,j,k));
#endif
                };
                auto bny = [&] (Array4<const Real> const& B) {
#if defined(WARPX_DIM_3D)
                    return 0.25_rt*(B(i-1,j,k-1) + B(i,j,k-1)
                                  + B(i-1,j,k) + B(i,j,k));
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    return 0.25_rt*(B(i-1,j-1,k) + B(i,j-1,k)
                                  + B(i-1,j,k) + B(i,j,k));
#else
                    return 0.5_rt*(B(i-1,j,k) + B(i,j,k));
#endif
                };
                auto bnz = [&] (Array4<const Real> const& B) {
#if defined(WARPX_DIM_3D)
                    return 0.25_rt*(B(i-1,j-1,k) + B(i,j-1,k)
                                  + B(i-1,j,k) + B(i,j,k));
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    return 0.5_rt*(B(i-1,j,k) + B(i,j,k));
#else
                    return B(i,j,k);
#endif
                };
                Real bx_v = stag ? bnx(Bx) : Bx(i,j,k);
                Real by_v = stag ? bny(By) : By(i,j,k);
                Real bz_v = stag ? bnz(Bz) : Bz(i,j,k);
                if (add_external) {
                    bx_v += stag ? bnx(Bxe) : Bxe(i,j,k);
                    by_v += stag ? bny(Bye) : Bye(i,j,k);
                    bz_v += stag ? bnz(Bze) : Bze(i,j,k);
                }
                const Real Bmag = std::sqrt(bx_v*bx_v + by_v*by_v + bz_v*bz_v);
                const Real rho_raw = r(i,j,k,0);
                const Real rho_e = HybridSmoothFloor(rho_raw, rho_floor, floor_w);
                if (Bmag > 0.0_rt) {
                    coef(i,j,k,0) = bx_v / Bmag;
                    coef(i,j,k,1) = by_v / Bmag;
                    coef(i,j,k,2) = bz_v / Bmag;
                } else {
                    coef(i,j,k,0) = 0.0_rt;
                    coef(i,j,k,1) = 0.0_rt;
                    coef(i,j,k,2) = 0.0_rt;
                }
                coef(i,j,k,3) = theta_dt * Bmag * one_over_mu0 / rho_e;
                // Resistivity frozen at the linearization state |J0|
                // (|J| = 0 would underestimate anomalous eta(J) in current
                // sheets). The parsers see the raw (unfloored) rho, matching
                // the Ohm's-law solve: flooring would deselect vacuum
                // branches if(rho < rho_vac, ...).
                const Real J0mag = std::sqrt(cJ(i,j,k,0)*cJ(i,j,k,0)
                                           + cJ(i,j,k,1)*cJ(i,j,k,1)
                                           + cJ(i,j,k,2)*cJ(i,j,k,2)) * one_over_mu0;
                coef(i,j,k,4) = theta_dt * eta(rho_raw, J0mag, t_now) * one_over_mu0;
                coef(i,j,k,5) = include_hyper ?
                    theta_dt * eta_h(rho_raw, Bmag) * one_over_mu0 : 0.0_rt;
                // ion-response strengths w_c = Sigma_cc |B_tot| / rho_c
                // (S v = (Sigma v) x B / rho_c through the Hall term; the
                // raw dJ/dE deposit already carries the theta*dt of the
                // implicit velocity integration, so w ~ Omega_ci theta*dt)
                Real w6 = 0.0_rt, w7 = 0.0_rt, w8 = 0.0_rt;
                if (have_mm) {
                    const Real bfac = Bmag / rho_e;
                    w6 = Sxx(i,j,k,cxx) * bfac;
                    w7 = Syy(i,j,k,cyy) * bfac;
                    w8 = Szz(i,j,k,czz) * bfac;
                }
                coef(i,j,k,6) = w6;
                coef(i,j,k,7) = w7;
                coef(i,j,k,8) = w8;
                coef(i,j,k,9) = rho_e;
                // pe row: pointwise diagonal D_p = 1 + coef(10), coef(10) =
                // theta*dt*[gamma*sum_d|ue_d|/dx_d + conduction]. |ue| from
                // the floored density: in floored cells the fictitious ue
                // only GROWS the diagonal, further damping the correction
                // there -- the safe direction for a PC.
                if (pe_row_on) {
                    const Real irho = 1.0_rt / rho_e;
                    // |ue_d| clamped like the residual's pe advance
                    auto uec = [=] (Real u) {
                        return (ue_cap_pc > 0.0_rt) ? ue_cap_pc * std::tanh(u / ue_cap_pc) : u;
                    };
                    Real adv = 0.0_rt;
#if defined(WARPX_DIM_1D_Z)
                    adv = uec(std::abs(Jiz(i,j,k) - Jpz(i,j,k)) * irho) * dxi_pc[0];
#elif defined(WARPX_DIM_XZ)
                    adv = uec(std::abs(Jix(i,j,k) - Jpx(i,j,k)) * irho) * dxi_pc[0]
                        + uec(std::abs(Jiz(i,j,k) - Jpz(i,j,k)) * irho) * dxi_pc[1];
#elif defined(WARPX_DIM_3D)
                    adv = uec(std::abs(Jix(i,j,k) - Jpx(i,j,k)) * irho) * dxi_pc[0]
                        + uec(std::abs(Jiy(i,j,k) - Jpy(i,j,k)) * irho) * dxi_pc[1]
                        + uec(std::abs(Jiz(i,j,k) - Jpz(i,j,k)) * irho) * dxi_pc[2];
#endif
                    // + the electron-conduction stiffness (kappa_e grad(pe/n)
                    // diffusion diagonal) so D_p tracks the pe row when
                    // conduction dominates the grid scale
                    coef(i,j,k,10) = theta_dt * (gamma_e * adv
                        + (gamma_e - 1.0_rt) * kappa_pc
                          * PhysConst::q_e * l2_pc / rho_e);
                } else {
                    coef(i,j,k,10) = 0.0_rt;
                }
                // pair: frozen electron velocity, stored pre-scaled as
                //   theta*dt * ue0 = theta*dt (J_i - J_plasma)/rho_c.
                // staggered: each J-like component is down-averaged to the
                // node over its cell-centered dimension (the residual's
                // ghost-filled J provides the lo-side reads), matching the
                // B down-averages above
                if (ue_on) {
                    const Real irho = theta_dt / rho_e;
                    if (stag) {
                        auto jn = [&] (Array4<const Real> const& Ja,
                                       Array4<const Real> const& Jb,
                                       int c) -> Real {
                            amrex::ignore_unused(c);
#if defined(WARPX_DIM_3D)
                            if (c == 0) {
                                return 0.5_rt*(Ja(i-1,j,k) - Jb(i-1,j,k)
                                             + Ja(i,j,k)   - Jb(i,j,k));
                            }
                            if (c == 1) {
                                return 0.5_rt*(Ja(i,j-1,k) - Jb(i,j-1,k)
                                             + Ja(i,j,k)   - Jb(i,j,k));
                            }
                            return 0.5_rt*(Ja(i,j,k-1) - Jb(i,j,k-1)
                                         + Ja(i,j,k)   - Jb(i,j,k));
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                            if (c == 0) {
                                return 0.5_rt*(Ja(i-1,j,k) - Jb(i-1,j,k)
                                             + Ja(i,j,k)   - Jb(i,j,k));
                            }
                            if (c == 2) {
                                return 0.5_rt*(Ja(i,j-1,k) - Jb(i,j-1,k)
                                             + Ja(i,j,k)   - Jb(i,j,k));
                            }
                            return Ja(i,j,k) - Jb(i,j,k);
#else
                            if (c == 2) {
                                return 0.5_rt*(Ja(i-1,j,k) - Jb(i-1,j,k)
                                             + Ja(i,j,k)   - Jb(i,j,k));
                            }
                            return Ja(i,j,k) - Jb(i,j,k);
#endif
                        };
                        coef(i,j,k,11) = jn(Jix, Jpx, 0) * irho;
                        coef(i,j,k,12) = jn(Jiy, Jpy, 1) * irho;
                        coef(i,j,k,13) = jn(Jiz, Jpz, 2) * irho;
                    } else {
                        coef(i,j,k,11) = (Jix(i,j,k) - Jpx(i,j,k)) * irho;
                        coef(i,j,k,12) = (Jiy(i,j,k) - Jpy(i,j,k)) * irho;
                        coef(i,j,k,13) = (Jiz(i,j,k) - Jpz(i,j,k)) * irho;
                    }
                } else {
                    coef(i,j,k,11) = 0.0_rt;
                    coef(i,j,k,12) = 0.0_rt;
                    coef(i,j,k,13) = 0.0_rt;
                }
                // pair: LOCAL W-hat scale s_W = max(beta_loc, 1)/l2 =
                // max(c_w, 1/l2), so every node's W-residual is weighted by
                // its own E-row coupling ||U_loc|| (a single global scale
                // cannot weight both the vacuum, beta ~ 1/rho_floor, and the
                // plasma of a strong density-contrast profile)
                coef(i,j,k,14) = amrex::max(coef(i,j,k,3), Real(1.0)/Real(l2_pc));
            });
    }

    // Peak whistler CFL: calibrates the inner-GMRES W-space tolerance so
    // it targets E-space accuracy (r_E = -U r_W; see the class docs)
    if (m_nested) {
        m_beta_max = static_cast<RT>(m_l2[0] * m_coef[0]->max(3, 0));
        if (m_staggered) {
            // On the staggered grid the E-space amplification of U is
            // carried by the interpolation-damped Hall chain: each
            // staggered member of the active transverse pair contributes
            // cos(theta/2) ~ 1/2 at the grid scale, so ||U|| (and with it
            // the tolerance scaling) is reduced by the least-damped pair's
            // factor -- both members staggered in 3D, one nodal member in
            // 2D/RZ (the undamped 1D pair keeps the full beta).
#if defined(WARPX_DIM_3D)
            m_beta_max *= RT(0.25);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            m_beta_max *= RT(0.5);
#endif
        }
    }

    // Coarse operators by rediscretization: inject the coefficients down the
    // hierarchy (l2 is recomputed per level in Define()) and exchange their
    // internal ghosts for the U-stage.
    FillGhosts(*m_coef[0], 0);
    ApplyDomainClosure(*m_coef[0], 0, ClosureClass::Coef);
    for (int lev = 1; lev < m_nlevels; lev++) {
        RestrictCoef(lev);
        FillGhosts(*m_coef[lev], lev);
        ApplyDomainClosure(*m_coef[lev], lev, ClosureClass::Coef);
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::RestrictCoef ( int a_lev )
{
    using namespace amrex;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_coef[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& c = m_coef[a_lev]->array(mfi);
        Array4<const Real> const& f = m_coef[a_lev-1]->const_array(mfi);
        amrex::ParallelFor(bx, kNCoef - 1,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int m)
            {
                // comps 0-9 (bhat, c_w, a_e, a_H, MM, rho) + the pair's
                // electron velocity (11-13) and local W-hat scale (14);
                // the level-0-only pe-row diagonal (10) is skipped
                const int n = (m < 10) ? m : (m + 1);
#if defined(WARPX_DIM_3D)
                c(i,j,k,n) = f(2*i,2*j,2*k,n);
#elif defined(WARPX_DIM_XZ)
                c(i,j,k,n) = f(2*i,2*j,k,n);
#else // WARPX_DIM_1D_Z
                c(i,j,k,n) = f(2*i,j,k,n);
#endif
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ComputeCurl ( const MF& a, MF& out, int a_lev, int a_grow,
                                       CurlKind a_kind, int a_scomp ) const
{
    using namespace amrex;

#if defined(WARPX_DIM_RZ)
    const auto* dx = m_geom[a_lev].CellSize();
    const Real dr = dx[0];
    const Real inv_dr = 1.0_rt/dx[0];
    const Real inv_dz = 1.0_rt/dx[1];
    const Real rmin = static_cast<Real>(m_geom[a_lev].ProbLo(0));
    const bool faraday = (a_kind == CurlKind::EtoB);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(a_grow);
        Array4<const Real> const& av = a.const_array(mfi, a_scomp);
        Array4<Real> const& cv = out.array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::ignore_unused(k);
                Real cr, ct, cz;
                if (faraday) {
                    hybridpcrz::CurlF(av, cr, ct, cz, i, j, rmin, dr,
                                      inv_dr, inv_dz);
                } else {
                    hybridpcrz::CurlA(av, cr, ct, cz, i, j, rmin, dr,
                                      inv_dr, inv_dz);
                }
                cv(i,j,0,0) = cr;
                cv(i,j,0,1) = ct;
                cv(i,j,0,2) = cz;
            });
    }
#else
    // staggered: compact up (E->B) / down (B->J) index-space differences
    const bool faraday = (a_kind == CurlKind::EtoB);
    const bool stag = m_staggered;
    const auto inv2dx = m_inv2dx[a_lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(a_grow);
        Array4<const Real> const& av = a.const_array(mfi, a_scomp);
        Array4<Real> const& cv = out.array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
#if defined(WARPX_DIM_3D)
                if (stag) {
                    const int  s  = faraday ? 0 : -1;
                    const Real ix = 2.0_rt*inv2dx[0];
                    const Real iy = 2.0_rt*inv2dx[1];
                    const Real iz = 2.0_rt*inv2dx[2];
                    cv(i,j,k,0) = (av(i,j+1+s,k,2) - av(i,j+s,k,2))*iy
                                - (av(i,j,k+1+s,1) - av(i,j,k+s,1))*iz;
                    cv(i,j,k,1) = (av(i,j,k+1+s,0) - av(i,j,k+s,0))*iz
                                - (av(i+1+s,j,k,2) - av(i+s,j,k,2))*ix;
                    cv(i,j,k,2) = (av(i+1+s,j,k,1) - av(i+s,j,k,1))*ix
                                - (av(i,j+1+s,k,0) - av(i,j+s,k,0))*iy;
                } else {
                    const Real dyaz = (av(i,j+1,k,2) - av(i,j-1,k,2)) * inv2dx[1];
                    const Real dzay = (av(i,j,k+1,1) - av(i,j,k-1,1)) * inv2dx[2];
                    const Real dzax = (av(i,j,k+1,0) - av(i,j,k-1,0)) * inv2dx[2];
                    const Real dxaz = (av(i+1,j,k,2) - av(i-1,j,k,2)) * inv2dx[0];
                    const Real dxay = (av(i+1,j,k,1) - av(i-1,j,k,1)) * inv2dx[0];
                    const Real dyax = (av(i,j+1,k,0) - av(i,j-1,k,0)) * inv2dx[1];
                    cv(i,j,k,0) = dyaz - dzay;
                    cv(i,j,k,1) = dzax - dxaz;
                    cv(i,j,k,2) = dxay - dyax;
                }
#elif defined(WARPX_DIM_XZ)
                // (i,j) = (x,z); d/dy = 0
                if (stag) {
                    const int  s  = faraday ? 0 : -1;
                    const Real ix = 2.0_rt*inv2dx[0];
                    const Real iz = 2.0_rt*inv2dx[2];
                    cv(i,j,k,0) = -(av(i,j+1+s,k,1) - av(i,j+s,k,1))*iz;
                    cv(i,j,k,1) = (av(i,j+1+s,k,0) - av(i,j+s,k,0))*iz
                                - (av(i+1+s,j,k,2) - av(i+s,j,k,2))*ix;
                    cv(i,j,k,2) = (av(i+1+s,j,k,1) - av(i+s,j,k,1))*ix;
                } else {
                    const Real dzay = (av(i,j+1,k,1) - av(i,j-1,k,1)) * inv2dx[2];
                    const Real dzax = (av(i,j+1,k,0) - av(i,j-1,k,0)) * inv2dx[2];
                    const Real dxaz = (av(i+1,j,k,2) - av(i-1,j,k,2)) * inv2dx[0];
                    const Real dxay = (av(i+1,j,k,1) - av(i-1,j,k,1)) * inv2dx[0];
                    cv(i,j,k,0) = -dzay;
                    cv(i,j,k,1) = dzax - dxaz;
                    cv(i,j,k,2) = dxay;
                }
#else // WARPX_DIM_1D_Z
                // i = z; d/dx = d/dy = 0
                if (stag) {
                    const int  s  = faraday ? 0 : -1;
                    const Real iz = 2.0_rt*inv2dx[2];
                    cv(i,j,k,0) = -(av(i+1+s,j,k,1) - av(i+s,j,k,1))*iz;
                    cv(i,j,k,1) = (av(i+1+s,j,k,0) - av(i+s,j,k,0))*iz;
                    cv(i,j,k,2) = 0.0_rt;
                } else {
                    const Real dzay = (av(i+1,j,k,1) - av(i-1,j,k,1)) * inv2dx[2];
                    const Real dzax = (av(i+1,j,k,0) - av(i-1,j,k,0)) * inv2dx[2];
                    cv(i,j,k,0) = -dzay;
                    cv(i,j,k,1) = dzax;
                    cv(i,j,k,2) = 0.0_rt;
                }
#endif
            });
    }
#endif // WARPX_DIM_RZ
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyUValid ( const MF& a, MF& out, int a_lev,
                                       int a_grow ) const
{
    using namespace amrex;

    const auto invdx2 = m_invdx2[a_lev];
    const bool stag = m_staggered;
#if defined(WARPX_DIM_RZ)
    amrex::ignore_unused(invdx2);   // the RZ Laplacian is hybridpcrz::LapJ
    const auto* dxg = m_geom[a_lev].CellSize();
    const Real rz_dr = dxg[0];
    const Real rz_inv_dr = 1.0_rt/dxg[0];
    const Real rz_inv_dz = 1.0_rt/dxg[1];
    const Real rz_rmin = static_cast<Real>(m_geom[a_lev].ProbLo(0));
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(a_grow);
        Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
        Array4<const Real> const& av = a.const_array(mfi);
        Array4<Real> const& ov = out.array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                const Real bhx = coef(i,j,k,0);
                const Real bhy = coef(i,j,k,1);
                const Real bhz = coef(i,j,k,2);
                const Real cw = coef(i,j,k,3);
                const Real ae = coef(i,j,k,4);
                const Real ah = coef(i,j,k,5);

                const Real ax = av(i,j,k,0);
                const Real ay = av(i,j,k,1);
                const Real az = av(i,j,k,2);

                Real lx = 0.0_rt, ly = 0.0_rt, lz = 0.0_rt;
#if defined(WARPX_DIM_3D)
                lx = (av(i+1,j,k,0) - 2.0_rt*ax + av(i-1,j,k,0))*invdx2[0]
                   + (av(i,j+1,k,0) - 2.0_rt*ax + av(i,j-1,k,0))*invdx2[1]
                   + (av(i,j,k+1,0) - 2.0_rt*ax + av(i,j,k-1,0))*invdx2[2];
                ly = (av(i+1,j,k,1) - 2.0_rt*ay + av(i-1,j,k,1))*invdx2[0]
                   + (av(i,j+1,k,1) - 2.0_rt*ay + av(i,j-1,k,1))*invdx2[1]
                   + (av(i,j,k+1,1) - 2.0_rt*ay + av(i,j,k-1,1))*invdx2[2];
                lz = (av(i+1,j,k,2) - 2.0_rt*az + av(i-1,j,k,2))*invdx2[0]
                   + (av(i,j+1,k,2) - 2.0_rt*az + av(i,j-1,k,2))*invdx2[1]
                   + (av(i,j,k+1,2) - 2.0_rt*az + av(i,j,k-1,2))*invdx2[2];
#elif defined(WARPX_DIM_XZ)
                lx = (av(i+1,j,k,0) - 2.0_rt*ax + av(i-1,j,k,0))*invdx2[0]
                   + (av(i,j+1,k,0) - 2.0_rt*ax + av(i,j-1,k,0))*invdx2[2];
                ly = (av(i+1,j,k,1) - 2.0_rt*ay + av(i-1,j,k,1))*invdx2[0]
                   + (av(i,j+1,k,1) - 2.0_rt*ay + av(i,j-1,k,1))*invdx2[2];
                lz = (av(i+1,j,k,2) - 2.0_rt*az + av(i-1,j,k,2))*invdx2[0]
                   + (av(i,j+1,k,2) - 2.0_rt*az + av(i,j-1,k,2))*invdx2[2];
#elif defined(WARPX_DIM_RZ)
                hybridpcrz::LapJ(av, lx, ly, lz, i, j, rz_rmin, rz_dr,
                                 rz_inv_dr, rz_inv_dz);
#else
                lx = (av(i+1,j,k,0) - 2.0_rt*ax + av(i-1,j,k,0))*invdx2[2];
                ly = (av(i+1,j,k,1) - 2.0_rt*ay + av(i-1,j,k,1))*invdx2[2];
                lz = (av(i+1,j,k,2) - 2.0_rt*az + av(i-1,j,k,2))*invdx2[2];
#endif
                if (stag) {
                    Real ux, uy, uz;
                    HybridPICPCUStag(av, coef, i, j, k, lx, ly, lz,
                                     ux, uy, uz);
                    amrex::ignore_unused(bhx, bhy, bhz, cw, ae, ah);
                    ov(i,j,k,0) = ux;
                    ov(i,j,k,1) = uy;
                    ov(i,j,k,2) = uz;
                } else {
                    ov(i,j,k,0) = cw*(ay*bhz - az*bhy) + ae*ax - ah*lx;
                    ov(i,j,k,1) = cw*(az*bhx - ax*bhz) + ae*ay - ah*ly;
                    ov(i,j,k,2) = cw*(ax*bhy - ay*bhx) + ae*az - ah*lz;
                }
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::StageKU ( const MF& a, int a_lev )
{
    using namespace amrex;

    const auto invdx2 = m_invdx2[a_lev];
    const bool ion = m_have_mm;
    const bool stag = m_staggered;
#if defined(WARPX_DIM_RZ)
    amrex::ignore_unused(invdx2);   // the RZ Laplacian is hybridpcrz::LapJ
    const auto* dxg = m_geom[a_lev].CellSize();
    const Real rz_dr = dxg[0];
    const Real rz_inv_dr = 1.0_rt/dxg[0];
    const Real rz_inv_dz = 1.0_rt/dxg[1];
    const Real rz_rmin = static_cast<Real>(m_geom[a_lev].ProbLo(0));
#endif

    // tmp = U a on growntilebox(kNGrow - 1): U reads a's first ghost ring
    // beyond that, so a must carry kNGrow filled ghosts
    const int gU = kNGrow - 1;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_tmp[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(gU);
        Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
        Array4<const Real> const& av = a.const_array(mfi);
        Array4<Real> const& ov = m_tmp[a_lev]->array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                const Real bhx = coef(i,j,k,0);
                const Real bhy = coef(i,j,k,1);
                const Real bhz = coef(i,j,k,2);
                const Real cw = coef(i,j,k,3);
                const Real ae = coef(i,j,k,4);
                const Real ah = coef(i,j,k,5);

                const Real ax = av(i,j,k,0);
                const Real ay = av(i,j,k,1);
                const Real az = av(i,j,k,2);

                Real lx = 0.0_rt, ly = 0.0_rt, lz = 0.0_rt;
#if defined(WARPX_DIM_3D)
                lx = (av(i+1,j,k,0) - 2.0_rt*ax + av(i-1,j,k,0))*invdx2[0]
                   + (av(i,j+1,k,0) - 2.0_rt*ax + av(i,j-1,k,0))*invdx2[1]
                   + (av(i,j,k+1,0) - 2.0_rt*ax + av(i,j,k-1,0))*invdx2[2];
                ly = (av(i+1,j,k,1) - 2.0_rt*ay + av(i-1,j,k,1))*invdx2[0]
                   + (av(i,j+1,k,1) - 2.0_rt*ay + av(i,j-1,k,1))*invdx2[1]
                   + (av(i,j,k+1,1) - 2.0_rt*ay + av(i,j,k-1,1))*invdx2[2];
                lz = (av(i+1,j,k,2) - 2.0_rt*az + av(i-1,j,k,2))*invdx2[0]
                   + (av(i,j+1,k,2) - 2.0_rt*az + av(i,j-1,k,2))*invdx2[1]
                   + (av(i,j,k+1,2) - 2.0_rt*az + av(i,j,k-1,2))*invdx2[2];
#elif defined(WARPX_DIM_XZ)
                lx = (av(i+1,j,k,0) - 2.0_rt*ax + av(i-1,j,k,0))*invdx2[0]
                   + (av(i,j+1,k,0) - 2.0_rt*ax + av(i,j-1,k,0))*invdx2[2];
                ly = (av(i+1,j,k,1) - 2.0_rt*ay + av(i-1,j,k,1))*invdx2[0]
                   + (av(i,j+1,k,1) - 2.0_rt*ay + av(i,j-1,k,1))*invdx2[2];
                lz = (av(i+1,j,k,2) - 2.0_rt*az + av(i-1,j,k,2))*invdx2[0]
                   + (av(i,j+1,k,2) - 2.0_rt*az + av(i,j-1,k,2))*invdx2[2];
#elif defined(WARPX_DIM_RZ)
                hybridpcrz::LapJ(av, lx, ly, lz, i, j, rz_rmin, rz_dr,
                                 rz_inv_dr, rz_inv_dz);
#else
                lx = (av(i+1,j,k,0) - 2.0_rt*ax + av(i-1,j,k,0))*invdx2[2];
                ly = (av(i+1,j,k,1) - 2.0_rt*ay + av(i-1,j,k,1))*invdx2[2];
                lz = (av(i+1,j,k,2) - 2.0_rt*az + av(i-1,j,k,2))*invdx2[2];
#endif
                Real ux, uy, uz;
                if (stag) {
                    // Yee: Hall rotation through the residual's node-
                    // interpolation chain, resistive legs at the component
                    // locations (see HybridPICPCUStag)
                    HybridPICPCUStag(av, coef, i, j, k, lx, ly, lz,
                                     ux, uy, uz);
                    amrex::ignore_unused(bhx, bhy, bhz, cw, ae, ah);
                } else {
                    ux = cw*(ay*bhz - az*bhy) + ae*ax - ah*lx;
                    uy = cw*(az*bhx - ax*bhz) + ae*ay - ah*ly;
                    uz = cw*(ax*bhy - ay*bhx) + ae*az - ah*lz;
                }
                // the W-system carries (I+S)^-1 U (deposited-mass-matrix
                // ion response, see class docs); identity when w = 0
                if (ion) {
                    HybridPICPCInvIpS(bhx, bhy, bhz,
                                      coef(i,j,k,6), coef(i,j,k,7),
                                      coef(i,j,k,8), ux, uy, uz);
                }
                ov(i,j,k,0) = ux;
                ov(i,j,k,1) = uy;
                ov(i,j,k,2) = uz;
            });
    }

    // The residual BCs the linearized Ohm output before the Faraday curl
    // (ApplyEfieldBoundary on E_ohm) and the resulting delta-B before the
    // Ampere curl (ApplyBfieldBoundary after EvolveB): mirror both here.
    ApplyDomainClosure(*m_tmp[a_lev], a_lev, ClosureClass::EField);

    // tmp2 = curl tmp on growntilebox(gU - 1)
    ComputeCurl(*m_tmp[a_lev], *m_tmp2[a_lev], a_lev, gU - 1, CurlKind::EtoB);
    ApplyDomainClosure(*m_tmp2[a_lev], a_lev, ClosureClass::BField);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::StagePair ( int a_lev )
{
    // F = curl xE on grown(gU), so consumers can read its first ring for
    // the inline kk = sW curl F
    const int gU = kNGrow - 1;
    ComputeCurl(*m_xF[a_lev], *m_tmp[a_lev], a_lev, gU, CurlKind::EtoB);
    ApplyDomainClosure(*m_tmp[a_lev], a_lev, ClosureClass::BField);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::SmoothLevel ( int a_lev, int a_nsweeps, bool a_x_is_zero )
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::SmoothLevel()");
#if defined(WARPX_DIM_RZ)
    const auto* dxg = m_geom[a_lev].CellSize();
    const Real rz_dr = dxg[0];
    const Real rz_inv_dr = 1.0_rt/dxg[0];
    const Real rz_inv_dz = 1.0_rt/dxg[1];
    const Real rz_rmin = static_cast<Real>(m_geom[a_lev].ProbLo(0));
#endif

    const RT l2 = m_l2[a_lev];
    const RT sigma = m_sigma;
    const bool stag = m_staggered;
    const auto inv2dx = m_inv2dx[a_lev];
#if defined(WARPX_DIM_RZ)
    amrex::ignore_unused(inv2dx);   // the RZ curl is hybridpcrz::CurlA
#endif

    for (int sweep = 0; sweep < a_nsweeps; sweep++) {

        const bool r_is_b = (sweep == 0 && a_x_is_zero);
        if (!r_is_b) {
            // one exchange feeds the whole U -> curl -> curl composition
            FillGhosts(*m_xW[a_lev], a_lev);
            // W is delta-J-like: parity ghosts, no wall write
            ApplyDomainClosure(*m_xW[a_lev], a_lev, ClosureClass::JField);
            if (m_pair) {
                // E-half of the primal pair
                FillGhosts(*m_xF[a_lev], a_lev);
                ApplyDomainClosure(*m_xF[a_lev], a_lev, ClosureClass::EField);
                StagePair(a_lev);
            } else {
                StageKU(*m_xW[a_lev], a_lev);
            }
        }

        if (m_pair) {
            // Joint per-node 2x2 block update (see the class docs): tc and
            // kk are evaluated inline (HybridPICPCPairTcKk) from the staged
            // F = curl xE and the xW neighbors.
            //   rE = bF - xE - tc, rW = bW - xW + kk
            //   P dE = rE - (1/sW) U_loc rW ; dW = rW + sW l2 dE
            const RT l1p = m_l1[a_lev];
            const auto invdx2p = m_invdx2[a_lev];
            const auto inv2dxp = m_inv2dx[a_lev];
            const bool ionp = m_have_mm;
            const Real wdl2 = static_cast<Real>(m_pair_wdf*l2);
            const RT sigmaw = m_sigma_w;
            if (stag) {
                // NODE-SANDWICH block: on the Yee grid the pointwise block's
                // cross-couplings would act between samples h/2 apart -- in
                // Fourier space a phase e^{i(theta_a-theta_b)/2} the real
                // rotation cannot represent. Every component-coupling term
                // of the update therefore goes through the operator's own
                // down-average/cross/up-average chain (symbol c_a c_b: real,
                // damped, phase-free); every diagonal term stays pointwise.
                //
                // stage A: exact residuals rE -> m_rF, rW -> m_rW (valid),
                //          then one exchange each fills the gather ghosts
                for (MFIter mfi(*m_coef[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
                    Array4<const Real> const& bW = m_bW[a_lev]->const_array(mfi);
                    Array4<const Real> const& bF = m_bF[a_lev]->const_array(mfi);
                    Array4<const Real> const& Fc = m_tmp[a_lev]->const_array(mfi);
                    Array4<const Real> const& xW = m_xW[a_lev]->const_array(mfi);
                    Array4<const Real> const& xE = m_xF[a_lev]->const_array(mfi);
                    Array4<Real> const& rE = m_rF[a_lev]->array(mfi);
                    Array4<Real> const& rW = m_rW[a_lev]->array(mfi);
                    amrex::ParallelFor(bx,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            if (r_is_b) {
                                for (int c = 0; c < 3; c++) {
                                    rE(i,j,k,c) = bF(i,j,k,c);
                                    rW(i,j,k,c) = bW(i,j,k,c);
                                }
                                return;
                            }
                            Real tcx, tcy, tcz, kkx, kky, kkz;
                            HybridPICPCPairTcKkStag(coef, xW, Fc, i, j, k, ionp,
                                invdx2p, inv2dxp,
                                tcx, tcy, tcz, kkx, kky, kkz);
                            rE(i,j,k,0) = bF(i,j,k,0) - xE(i,j,k,0) - tcx;
                            rE(i,j,k,1) = bF(i,j,k,1) - xE(i,j,k,1) - tcy;
                            rE(i,j,k,2) = bF(i,j,k,2) - xE(i,j,k,2) - tcz;
                            rW(i,j,k,0) = bW(i,j,k,0) - xW(i,j,k,0) + kkx;
                            rW(i,j,k,1) = bW(i,j,k,1) - xW(i,j,k,1) + kky;
                            rW(i,j,k,2) = bW(i,j,k,2) - xW(i,j,k,2) + kkz;
                        });
                }
                FillGhosts(*m_rF[a_lev], a_lev);
                FillGhosts(*m_rW[a_lev], a_lev);
                // stage B: node coupling q = Q sD~ + P^-1 sX~ -> m_tmp2 on
                // grow(1) (stage C reads q at the two nodes of each edge)
                for (MFIter mfi(*m_coef[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.growntilebox(1);
                    Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
                    Array4<const Real> const& rE = m_rF[a_lev]->const_array(mfi);
                    Array4<const Real> const& rW = m_rW[a_lev]->const_array(mfi);
                    Array4<Real> const& qq = m_tmp2[a_lev]->array(mfi);
                    amrex::ParallelFor(bx,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            // A- gather of both residual halves to the node
                            // (each component averaged along its own
                            // cell-centered dimension)
                            Real gEx, gEy, gEz, gWx, gWy, gWz;
#if defined(WARPX_DIM_3D)
                            gEx = 0.5_rt*(rE(i-1,j,k,0) + rE(i,j,k,0));
                            gEy = 0.5_rt*(rE(i,j-1,k,1) + rE(i,j,k,1));
                            gEz = 0.5_rt*(rE(i,j,k-1,2) + rE(i,j,k,2));
                            gWx = 0.5_rt*(rW(i-1,j,k,0) + rW(i,j,k,0));
                            gWy = 0.5_rt*(rW(i,j-1,k,1) + rW(i,j,k,1));
                            gWz = 0.5_rt*(rW(i,j,k-1,2) + rW(i,j,k,2));
#elif defined(WARPX_DIM_XZ)
                            gEx = 0.5_rt*(rE(i-1,j,k,0) + rE(i,j,k,0));
                            gEy = rE(i,j,k,1);
                            gEz = 0.5_rt*(rE(i,j-1,k,2) + rE(i,j,k,2));
                            gWx = 0.5_rt*(rW(i-1,j,k,0) + rW(i,j,k,0));
                            gWy = rW(i,j,k,1);
                            gWz = 0.5_rt*(rW(i,j-1,k,2) + rW(i,j,k,2));
#else
                            gEx = rE(i,j,k,0);
                            gEy = rE(i,j,k,1);
                            gEz = 0.5_rt*(rE(i-1,j,k,2) + rE(i,j,k,2));
                            gWx = rW(i,j,k,0);
                            gWy = rW(i,j,k,1);
                            gWz = 0.5_rt*(rW(i-1,j,k,2) + rW(i,j,k,2));
#endif
                            const Real bhx = coef(i,j,k,0);
                            const Real bhy = coef(i,j,k,1);
                            const Real bhz = coef(i,j,k,2);
                            const Real cw = coef(i,j,k,3);
                            const Real aU = coef(i,j,k,4) + coef(i,j,k,5)*l2;
                            const Real uem = std::sqrt(
                                coef(i,j,k,11)*coef(i,j,k,11)
                              + coef(i,j,k,12)*coef(i,j,k,12)
                              + coef(i,j,k,13)*coef(i,j,k,13));
                            const Real alpha = 1.0_rt + l1p*uem + l2*aU;
                            const Real beta = l2*cw;
                            const Real isW = 1.0_rt/coef(i,j,k,14);
                            // sD = gE - (1/sW) aU gW ; sX = -(1/sW) cw gW x b
                            const Real sDx = gEx - isW*aU*gWx;
                            const Real sDy = gEy - isW*aU*gWy;
                            const Real sDz = gEz - isW*aU*gWz;
                            const Real sXx = -isW*cw*(gWy*bhz - gWz*bhy);
                            const Real sXy = -isW*cw*(gWz*bhx - gWx*bhz);
                            const Real sXz = -isW*cw*(gWx*bhy - gWy*bhx);
                            // q = Q sD + P^-1 sX, with t = sD + sX:
                            //   q = [alpha sX - beta (t x b)
                            //        + (beta^2/alpha)(t.b) b
                            //        - (beta^2/alpha) sD] / (alpha^2+beta^2)
                            const Real tx = sDx + sXx;
                            const Real ty = sDy + sXy;
                            const Real tz = sDz + sXz;
                            const Real txb_x = ty*bhz - tz*bhy;
                            const Real txb_y = tz*bhx - tx*bhz;
                            const Real txb_z = tx*bhy - ty*bhx;
                            const Real tdb = tx*bhx + ty*bhy + tz*bhz;
                            const Real b2oa = beta*beta/alpha;
                            const Real inv_a2b2 = 1.0_rt/(alpha*alpha + beta*beta);
                            qq(i,j,k,0) = inv_a2b2*(alpha*sXx - beta*txb_x
                                          + b2oa*(tdb*bhx - sDx));
                            qq(i,j,k,1) = inv_a2b2*(alpha*sXy - beta*txb_y
                                          + b2oa*(tdb*bhy - sDy));
                            qq(i,j,k,2) = inv_a2b2*(alpha*sXz - beta*txb_z
                                          + b2oa*(tdb*bhz - sDz));
                        });
                }
                // stage C: dE = sD/alpha + A+ q ; xE += sigma dE,
                // xW += sigma (rW + sW l2 dE)
                for (MFIter mfi(*m_coef[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
                    Array4<const Real> const& rE = m_rF[a_lev]->const_array(mfi);
                    Array4<const Real> const& rW = m_rW[a_lev]->const_array(mfi);
                    Array4<const Real> const& qq = m_tmp2[a_lev]->const_array(mfi);
                    Array4<Real> const& xW = m_xW[a_lev]->array(mfi);
                    Array4<Real> const& xE = m_xF[a_lev]->array(mfi);
                    amrex::ParallelFor(bx,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            const Real aU = coef(i,j,k,4) + coef(i,j,k,5)*l2;
                            const Real uem = std::sqrt(
                                coef(i,j,k,11)*coef(i,j,k,11)
                              + coef(i,j,k,12)*coef(i,j,k,12)
                              + coef(i,j,k,13)*coef(i,j,k,13));
                            const Real ialpha = 1.0_rt/(1.0_rt + l1p*uem + l2*aU);
                            const Real swl = coef(i,j,k,14);
                            const Real isW = 1.0_rt/swl;
                            const Real l2s = static_cast<Real>(l2)*swl;
                            // A+ scatter of q to each component's edge
                            Real qx, qy, qz;
#if defined(WARPX_DIM_3D)
                            qx = 0.5_rt*(qq(i,j,k,0) + qq(i+1,j,k,0));
                            qy = 0.5_rt*(qq(i,j,k,1) + qq(i,j+1,k,1));
                            qz = 0.5_rt*(qq(i,j,k,2) + qq(i,j,k+1,2));
#elif defined(WARPX_DIM_XZ)
                            qx = 0.5_rt*(qq(i,j,k,0) + qq(i+1,j,k,0));
                            qy = qq(i,j,k,1);
                            qz = 0.5_rt*(qq(i,j,k,2) + qq(i,j+1,k,2));
#else
                            qx = qq(i,j,k,0);
                            qy = qq(i,j,k,1);
                            qz = 0.5_rt*(qq(i,j,k,2) + qq(i+1,j,k,2));
#endif
                            const Real dEx = (rE(i,j,k,0) - isW*aU*rW(i,j,k,0))*ialpha + qx;
                            const Real dEy = (rE(i,j,k,1) - isW*aU*rW(i,j,k,1))*ialpha + qy;
                            const Real dEz = (rE(i,j,k,2) - isW*aU*rW(i,j,k,2))*ialpha + qz;
                            xE(i,j,k,0) += sigma*dEx;
                            xE(i,j,k,1) += sigma*dEy;
                            xE(i,j,k,2) += sigma*dEz;
                            xW(i,j,k,0) += sigma*(rW(i,j,k,0) + l2s*dEx);
                            xW(i,j,k,1) += sigma*(rW(i,j,k,1) + l2s*dEy);
                            xW(i,j,k,2) += sigma*(rW(i,j,k,2) + l2s*dEz);
                        });
                }
                continue;
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*m_coef[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
                Array4<const Real> const& bW = m_bW[a_lev]->const_array(mfi);
                Array4<const Real> const& bF = m_bF[a_lev]->const_array(mfi);
                Array4<const Real> const& Fc = m_tmp[a_lev]->const_array(mfi);
                // reads from the current iterate, writes into the ping-pong
                // partner (a pure Jacobi sweep; see m_xW2)
                Array4<const Real> const& xW = m_xW[a_lev]->const_array(mfi);
                Array4<const Real> const& xE = m_xF[a_lev]->const_array(mfi);
                Array4<Real> const& xWo = m_xW2[a_lev]->array(mfi);
                Array4<Real> const& xEo = m_xF2[a_lev]->array(mfi);
                amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        Real rEx, rEy, rEz, rWx, rWy, rWz;
                        // current iterate; zero on the first sweep of a zero
                        // initial guess (the level's x is NOT pre-zeroed, see
                        // m_zero_free)
                        Real xEx = 0.0_rt, xEy = 0.0_rt, xEz = 0.0_rt;
                        Real xWx = 0.0_rt, xWy = 0.0_rt, xWz = 0.0_rt;
                        if (r_is_b) {
                            rEx = bF(i,j,k,0); rEy = bF(i,j,k,1); rEz = bF(i,j,k,2);
                            rWx = bW(i,j,k,0); rWy = bW(i,j,k,1); rWz = bW(i,j,k,2);
                        } else {
                            Real tcx, tcy, tcz, kkx, kky, kkz;
                            HybridPICPCPairTcKk(coef, xW, Fc, i, j, k, ionp,
                                invdx2p, inv2dxp,
                                tcx, tcy, tcz, kkx, kky, kkz);
                            xEx = xE(i,j,k,0); xEy = xE(i,j,k,1); xEz = xE(i,j,k,2);
                            xWx = xW(i,j,k,0); xWy = xW(i,j,k,1); xWz = xW(i,j,k,2);
                            rEx = bF(i,j,k,0) - xEx - tcx;
                            rEy = bF(i,j,k,1) - xEy - tcy;
                            rEz = bF(i,j,k,2) - xEz - tcz;
                            rWx = bW(i,j,k,0) - xWx + kkx;
                            rWy = bW(i,j,k,1) - xWy + kky;
                            rWz = bW(i,j,k,2) - xWz + kkz;
                        }
                        const Real bhx = coef(i,j,k,0);
                        const Real bhy = coef(i,j,k,1);
                        const Real bhz = coef(i,j,k,2);
                        const Real cw = coef(i,j,k,3);
                        const Real aU = coef(i,j,k,4) + coef(i,j,k,5)*l2;
                        // local W-hat scale (coef comp 14)
                        const Real swl = coef(i,j,k,14);
                        const Real isW = 1.0_rt/swl;
                        // s = rE - (1/sW) U_loc rW (the E-row couples to
                        // the scaled W-hat through U_loc/sW)
                        const Real sx = rEx - isW*(cw*(rWy*bhz - rWz*bhy) + aU*rWx);
                        const Real sy = rEy - isW*(cw*(rWz*bhx - rWx*bhz) + aU*rWy);
                        const Real sz = rEz - isW*(cw*(rWx*bhy - rWy*bhx) + aU*rWz);
                        // P^-1 s with the block defects: hyper-resistive at
                        // the Laplacian Nyquist symbol 2*l2, whistler at wdl2
                        const Real uem = std::sqrt(
                            coef(i,j,k,11)*coef(i,j,k,11)
                          + coef(i,j,k,12)*coef(i,j,k,12)
                          + coef(i,j,k,13)*coef(i,j,k,13));
                        const Real aUb = coef(i,j,k,4) + coef(i,j,k,5)*(2.0_rt*l2);
                        const Real alpha = 1.0_rt + l1p*uem + l2*aUb;
                        const Real beta = wdl2*cw;
                        const Real sxb_x = sy*bhz - sz*bhy;
                        const Real sxb_y = sz*bhx - sx*bhz;
                        const Real sxb_z = sx*bhy - sy*bhx;
                        const Real sdb = sx*bhx + sy*bhy + sz*bhz;
                        const Real inv_a2b2 = 1.0_rt / (alpha*alpha + beta*beta);
                        const Real b2oa = beta*beta/alpha;
                        const Real dEx = inv_a2b2*(alpha*sx - beta*sxb_x + b2oa*sdb*bhx);
                        const Real dEy = inv_a2b2*(alpha*sy - beta*sxb_y + b2oa*sdb*bhy);
                        const Real dEz = inv_a2b2*(alpha*sz - beta*sxb_z + b2oa*sdb*bhz);
                        xEo(i,j,k,0) = xEx + sigma*dEx;
                        xEo(i,j,k,1) = xEy + sigma*dEy;
                        xEo(i,j,k,2) = xEz + sigma*dEz;
                        // perp-b part of dE couples through the whistler
                        // defect, the par-b part through the full l2
                        const Real dEb = dEx*bhx + dEy*bhy + dEz*bhz;
                        const Real dl = static_cast<Real>(l2) - wdl2;
                        xWo(i,j,k,0) = xWx + sigmaw*(rWx + swl*(wdl2*dEx + dl*dEb*bhx));
                        xWo(i,j,k,1) = xWy + sigmaw*(rWy + swl*(wdl2*dEy + dl*dEb*bhy));
                        xWo(i,j,k,2) = xWz + sigmaw*(rWz + swl*(wdl2*dEz + dl*dEb*bhz));
                    });
            }
            // the partner now holds the new iterate (valid region); its
            // ghosts are refilled by the next FillGhosts before any read
            std::swap(m_xW[a_lev], m_xW2[a_lev]);
            std::swap(m_xF[a_lev], m_xF2[a_lev]);
        } else {
            // Staged residual + damped block update (in place):
            //   r = bW - xW - curl(tmp2)  (r = bW on the first zero-x sweep)
            //   xW += sigma P^-1 r,  P = alpha I + beta (. x bhat)
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*m_coef[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
                Array4<const Real> const& bW = m_bW[a_lev]->const_array(mfi);
                Array4<const Real> const& t2 = m_tmp2[a_lev]->const_array(mfi);
                Array4<Real> const& xW = m_xW[a_lev]->array(mfi);
                amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        Real sx, sy, sz;
                        if (r_is_b) {
                            sx = bW(i,j,k,0);
                            sy = bW(i,j,k,1);
                            sz = bW(i,j,k,2);
                        } else {
                            // curl(tmp2) inline
#if defined(WARPX_DIM_RZ)
                            Real ccx, ccy, ccz;
                            hybridpcrz::CurlA(t2, ccx, ccy, ccz, i, j, rz_rmin,
                                              rz_dr, rz_inv_dr, rz_inv_dz);
#elif defined(WARPX_DIM_3D)
                            Real ccx, ccy, ccz;
                            if (stag) {
                                const Real ix = 2.0_rt*inv2dx[0];
                                const Real iy = 2.0_rt*inv2dx[1];
                                const Real iz = 2.0_rt*inv2dx[2];
                                ccx = (t2(i,j,k,2) - t2(i,j-1,k,2))*iy
                                    - (t2(i,j,k,1) - t2(i,j,k-1,1))*iz;
                                ccy = (t2(i,j,k,0) - t2(i,j,k-1,0))*iz
                                    - (t2(i,j,k,2) - t2(i-1,j,k,2))*ix;
                                ccz = (t2(i,j,k,1) - t2(i-1,j,k,1))*ix
                                    - (t2(i,j,k,0) - t2(i,j-1,k,0))*iy;
                            } else {
                                ccx = (t2(i,j+1,k,2) - t2(i,j-1,k,2))*inv2dx[1]
                                    - (t2(i,j,k+1,1) - t2(i,j,k-1,1))*inv2dx[2];
                                ccy = (t2(i,j,k+1,0) - t2(i,j,k-1,0))*inv2dx[2]
                                    - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                                ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0]
                                    - (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[1];
                            }
#elif defined(WARPX_DIM_XZ)
                            Real ccx, ccy, ccz;
                            if (stag) {
                                const Real ix = 2.0_rt*inv2dx[0];
                                const Real iz = 2.0_rt*inv2dx[2];
                                ccx = -(t2(i,j,k,1) - t2(i,j-1,k,1))*iz;
                                ccy = (t2(i,j,k,0) - t2(i,j-1,k,0))*iz
                                    - (t2(i,j,k,2) - t2(i-1,j,k,2))*ix;
                                ccz = (t2(i,j,k,1) - t2(i-1,j,k,1))*ix;
                            } else {
                                ccx = -(t2(i,j+1,k,1) - t2(i,j-1,k,1))*inv2dx[2];
                                ccy = (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[2]
                                    - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                                ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0];
                            }
#else
                            Real ccx, ccy, ccz;
                            if (stag) {
                                const Real iz = 2.0_rt*inv2dx[2];
                                ccx = -(t2(i,j,k,1) - t2(i-1,j,k,1))*iz;
                                ccy = (t2(i,j,k,0) - t2(i-1,j,k,0))*iz;
                                ccz = 0.0_rt;
                            } else {
                                ccx = -(t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[2];
                                ccy = (t2(i+1,j,k,0) - t2(i-1,j,k,0))*inv2dx[2];
                                ccz = 0.0_rt;
                            }
#endif
                            sx = bW(i,j,k,0) - xW(i,j,k,0) - ccx;
                            sy = bW(i,j,k,1) - xW(i,j,k,1) - ccy;
                            sz = bW(i,j,k,2) - xW(i,j,k,2) - ccz;
                        }

                        const Real bhx = coef(i,j,k,0);
                        const Real bhy = coef(i,j,k,1);
                        const Real bhz = coef(i,j,k,2);
                        const Real cw = coef(i,j,k,3);
                        const Real aU = coef(i,j,k,4) + coef(i,j,k,5)*l2;

                        const Real alpha = 1.0_rt + l2*aU;
                        // ion-response contribution folded in through its
                        // isotropic (rotation) part, |S| ~ Omega_ci theta*dt.
                        // The full l2 c_w UPPER-bounds the rotation symbol
                        // on both grids (on the staggered grid the node-
                        // interpolation cosines damp the Hall chain pair-
                        // anisotropically, e.g. not at all for a nodal
                        // transverse pair): over-estimating beta keeps the
                        // damped sweep robust, at the cost of slow rotation
                        // modes the nested mode's Krylov space absorbs
                        const Real beta = l2*cw
                            + (coef(i,j,k,6) + coef(i,j,k,7) + coef(i,j,k,8))
                              * (1.0_rt/3.0_rt);
                        const Real inv_a2b2 = 1.0_rt / (alpha*alpha + beta*beta);
                        const Real b2oa = beta*beta/alpha;

                        const Real sxb_x = sy*bhz - sz*bhy;
                        const Real sxb_y = sz*bhx - sx*bhz;
                        const Real sxb_z = sx*bhy - sy*bhx;
                        const Real sdb = sx*bhx + sy*bhy + sz*bhz;
                        xW(i,j,k,0) += sigma*inv_a2b2*(alpha*sx - beta*sxb_x + b2oa*sdb*bhx);
                        xW(i,j,k,1) += sigma*inv_a2b2*(alpha*sy - beta*sxb_y + b2oa*sdb*bhy);
                        xW(i,j,k,2) += sigma*inv_a2b2*(alpha*sz - beta*sxb_z + b2oa*sdb*bhz);
                    });
            }
        }
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ResidualLevel ( int a_lev )
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::ResidualLevel()");
#if defined(WARPX_DIM_RZ)
    const auto* dxg = m_geom[a_lev].CellSize();
    const Real rz_dr = dxg[0];
    const Real rz_inv_dr = 1.0_rt/dxg[0];
    const Real rz_inv_dz = 1.0_rt/dxg[1];
    const Real rz_rmin = static_cast<Real>(m_geom[a_lev].ProbLo(0));
#endif

    const auto inv2dx = m_inv2dx[a_lev];
    const auto invdx2 = m_invdx2[a_lev];
    const bool stag = m_staggered;

    FillGhosts(*m_xW[a_lev], a_lev);
    ApplyDomainClosure(*m_xW[a_lev], a_lev, ClosureClass::JField);
    if (m_pair) {
        FillGhosts(*m_xF[a_lev], a_lev);
        ApplyDomainClosure(*m_xF[a_lev], a_lev, ClosureClass::EField);
        StagePair(a_lev);
        // pair residuals: rE = bF - xE - tc ; rW = bW - xW + kk
        const bool ionp = m_have_mm;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_rW[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> const& coef = m_coef[a_lev]->const_array(mfi);
            Array4<const Real> const& bW = m_bW[a_lev]->const_array(mfi);
            Array4<const Real> const& bF = m_bF[a_lev]->const_array(mfi);
            Array4<const Real> const& xW = m_xW[a_lev]->const_array(mfi);
            Array4<const Real> const& xE = m_xF[a_lev]->const_array(mfi);
            Array4<const Real> const& Fc = m_tmp[a_lev]->const_array(mfi);
            Array4<Real> const& rW = m_rW[a_lev]->array(mfi);
            Array4<Real> const& rF = m_rF[a_lev]->array(mfi);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real tcx, tcy, tcz, kkx, kky, kkz;
                    if (stag) {
                        HybridPICPCPairTcKkStag(coef, xW, Fc, i, j, k, ionp,
                                                invdx2, inv2dx,
                                                tcx, tcy, tcz, kkx, kky, kkz);
                    } else {
                        HybridPICPCPairTcKk(coef, xW, Fc, i, j, k, ionp,
                                            invdx2, inv2dx,
                                            tcx, tcy, tcz, kkx, kky, kkz);
                    }
                    rF(i,j,k,0) = bF(i,j,k,0) - xE(i,j,k,0) - tcx;
                    rF(i,j,k,1) = bF(i,j,k,1) - xE(i,j,k,1) - tcy;
                    rF(i,j,k,2) = bF(i,j,k,2) - xE(i,j,k,2) - tcz;
                    rW(i,j,k,0) = bW(i,j,k,0) - xW(i,j,k,0) + kkx;
                    rW(i,j,k,1) = bW(i,j,k,1) - xW(i,j,k,1) + kky;
                    rW(i,j,k,2) = bW(i,j,k,2) - xW(i,j,k,2) + kkz;
                });
        }
        return;
    }

    StageKU(*m_xW[a_lev], a_lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_rW[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<const Real> const& bW = m_bW[a_lev]->const_array(mfi);
        Array4<const Real> const& xW = m_xW[a_lev]->const_array(mfi);
        Array4<const Real> const& t2 = m_tmp2[a_lev]->const_array(mfi);
        Array4<Real> const& rW = m_rW[a_lev]->array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
#if defined(WARPX_DIM_RZ)
                Real ccx, ccy, ccz;
                hybridpcrz::CurlA(t2, ccx, ccy, ccz, i, j, rz_rmin,
                                  rz_dr, rz_inv_dr, rz_inv_dz);
#elif defined(WARPX_DIM_3D)
                Real ccx, ccy, ccz;
                if (stag) {
                    const Real ix = 2.0_rt*inv2dx[0];
                    const Real iy = 2.0_rt*inv2dx[1];
                    const Real iz = 2.0_rt*inv2dx[2];
                    ccx = (t2(i,j,k,2) - t2(i,j-1,k,2))*iy
                        - (t2(i,j,k,1) - t2(i,j,k-1,1))*iz;
                    ccy = (t2(i,j,k,0) - t2(i,j,k-1,0))*iz
                        - (t2(i,j,k,2) - t2(i-1,j,k,2))*ix;
                    ccz = (t2(i,j,k,1) - t2(i-1,j,k,1))*ix
                        - (t2(i,j,k,0) - t2(i,j-1,k,0))*iy;
                } else {
                    ccx = (t2(i,j+1,k,2) - t2(i,j-1,k,2))*inv2dx[1]
                        - (t2(i,j,k+1,1) - t2(i,j,k-1,1))*inv2dx[2];
                    ccy = (t2(i,j,k+1,0) - t2(i,j,k-1,0))*inv2dx[2]
                        - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                    ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0]
                        - (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[1];
                }
#elif defined(WARPX_DIM_XZ)
                Real ccx, ccy, ccz;
                if (stag) {
                    const Real ix = 2.0_rt*inv2dx[0];
                    const Real iz = 2.0_rt*inv2dx[2];
                    ccx = -(t2(i,j,k,1) - t2(i,j-1,k,1))*iz;
                    ccy = (t2(i,j,k,0) - t2(i,j-1,k,0))*iz
                        - (t2(i,j,k,2) - t2(i-1,j,k,2))*ix;
                    ccz = (t2(i,j,k,1) - t2(i-1,j,k,1))*ix;
                } else {
                    ccx = -(t2(i,j+1,k,1) - t2(i,j-1,k,1))*inv2dx[2];
                    ccy = (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[2]
                        - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                    ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0];
                }
#else
                Real ccx, ccy, ccz;
                if (stag) {
                    const Real iz = 2.0_rt*inv2dx[2];
                    ccx = -(t2(i,j,k,1) - t2(i-1,j,k,1))*iz;
                    ccy = (t2(i,j,k,0) - t2(i-1,j,k,0))*iz;
                    ccz = 0.0_rt;
                } else {
                    ccx = -(t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[2];
                    ccy = (t2(i+1,j,k,0) - t2(i-1,j,k,0))*inv2dx[2];
                    ccz = 0.0_rt;
                }
#endif
                rW(i,j,k,0) = bW(i,j,k,0) - xW(i,j,k,0) - ccx;
                rW(i,j,k,1) = bW(i,j,k,1) - xW(i,j,k,1) - ccy;
                rW(i,j,k,2) = bW(i,j,k,2) - xW(i,j,k,2) - ccz;
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyReducedOp ( const MF& a_blk, MF& a_out )
{
    using namespace amrex;
#if defined(WARPX_DIM_RZ)
    const auto* dxg = m_geom[0].CellSize();
    const Real rz_dr = dxg[0];
    const Real rz_inv_dr = 1.0_rt/dxg[0];
    const Real rz_inv_dz = 1.0_rt/dxg[1];
    const Real rz_rmin = static_cast<Real>(m_geom[0].ProbLo(0));
#endif

    // a_blk must be ghost-filled (kNGrow layers) by the caller; the W-half
    // is comps 0..2, the E-half (pair) comps 3..5
    const MF& a = a_blk;
    MF& out = a_out;
    if (m_pair) {
        // pair matvec: out_W = a - kk ; out_E = xE + tc, with F = curl xE
        // staged from the block's E-half
        const int gU = kNGrow - 1;
        ComputeCurl(a_blk, *m_tmp[0], 0, gU, CurlKind::EtoB, 3);
        ApplyDomainClosure(*m_tmp[0], 0, ClosureClass::BField);
        const bool ionp = m_have_mm;
        const bool stagp = m_staggered;
        const auto invdx2p = m_invdx2[0];
        const auto inv2dxp = m_inv2dx[0];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> const& coef = m_coef[0]->const_array(mfi);
            Array4<const Real> const& av = a.const_array(mfi);
            Array4<const Real> const& zE = a.const_array(mfi, 3);
            Array4<const Real> const& Fc = m_tmp[0]->const_array(mfi);
            Array4<Real> const& ow = out.array(mfi);
            Array4<Real> const& oe = out.array(mfi, 3);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real tcx, tcy, tcz, kkx, kky, kkz;
                    if (stagp) {
                        HybridPICPCPairTcKkStag(coef, av, Fc, i, j, k, ionp,
                                                invdx2p, inv2dxp,
                                                tcx, tcy, tcz, kkx, kky, kkz);
                    } else {
                        HybridPICPCPairTcKk(coef, av, Fc, i, j, k, ionp,
                                            invdx2p, inv2dxp,
                                            tcx, tcy, tcz, kkx, kky, kkz);
                    }
                    ow(i,j,k,0) = av(i,j,k,0) - kkx;
                    ow(i,j,k,1) = av(i,j,k,1) - kky;
                    ow(i,j,k,2) = av(i,j,k,2) - kkz;
                    oe(i,j,k,0) = zE(i,j,k,0) + tcx;
                    oe(i,j,k,1) = zE(i,j,k,1) + tcy;
                    oe(i,j,k,2) = zE(i,j,k,2) + tcz;
                });
        }
        return;
    }

    // reduced W-system: out = a + curl curl (U a)
    StageKU(a, 0);
    const auto inv2dx = m_inv2dx[0];
    const bool stag = m_staggered;
#if defined(WARPX_DIM_RZ)
    amrex::ignore_unused(inv2dx, stag);   // the RZ curl is hybridpcrz::CurlA
#endif
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<const Real> const& av = a.const_array(mfi);
        Array4<const Real> const& t2 = m_tmp2[0]->const_array(mfi);
        Array4<Real> const& ov = out.array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
#if defined(WARPX_DIM_RZ)
                Real ccx, ccy, ccz;
                hybridpcrz::CurlA(t2, ccx, ccy, ccz, i, j, rz_rmin,
                                  rz_dr, rz_inv_dr, rz_inv_dz);
#elif defined(WARPX_DIM_3D)
                Real ccx, ccy, ccz;
                if (stag) {
                    const Real ix = 2.0_rt*inv2dx[0];
                    const Real iy = 2.0_rt*inv2dx[1];
                    const Real iz = 2.0_rt*inv2dx[2];
                    ccx = (t2(i,j,k,2) - t2(i,j-1,k,2))*iy
                        - (t2(i,j,k,1) - t2(i,j,k-1,1))*iz;
                    ccy = (t2(i,j,k,0) - t2(i,j,k-1,0))*iz
                        - (t2(i,j,k,2) - t2(i-1,j,k,2))*ix;
                    ccz = (t2(i,j,k,1) - t2(i-1,j,k,1))*ix
                        - (t2(i,j,k,0) - t2(i,j-1,k,0))*iy;
                } else {
                    ccx = (t2(i,j+1,k,2) - t2(i,j-1,k,2))*inv2dx[1]
                        - (t2(i,j,k+1,1) - t2(i,j,k-1,1))*inv2dx[2];
                    ccy = (t2(i,j,k+1,0) - t2(i,j,k-1,0))*inv2dx[2]
                        - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                    ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0]
                        - (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[1];
                }
#elif defined(WARPX_DIM_XZ)
                Real ccx, ccy, ccz;
                if (stag) {
                    const Real ix = 2.0_rt*inv2dx[0];
                    const Real iz = 2.0_rt*inv2dx[2];
                    ccx = -(t2(i,j,k,1) - t2(i,j-1,k,1))*iz;
                    ccy = (t2(i,j,k,0) - t2(i,j-1,k,0))*iz
                        - (t2(i,j,k,2) - t2(i-1,j,k,2))*ix;
                    ccz = (t2(i,j,k,1) - t2(i-1,j,k,1))*ix;
                } else {
                    ccx = -(t2(i,j+1,k,1) - t2(i,j-1,k,1))*inv2dx[2];
                    ccy = (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[2]
                        - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                    ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0];
                }
#else
                Real ccx, ccy, ccz;
                if (stag) {
                    const Real iz = 2.0_rt*inv2dx[2];
                    ccx = -(t2(i,j,k,1) - t2(i-1,j,k,1))*iz;
                    ccy = (t2(i,j,k,0) - t2(i-1,j,k,0))*iz;
                    ccz = 0.0_rt;
                } else {
                    ccx = -(t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[2];
                    ccy = (t2(i+1,j,k,0) - t2(i-1,j,k,0))*inv2dx[2];
                    ccz = 0.0_rt;
                }
#endif
                ov(i,j,k,0) = av(i,j,k,0) + ccx;
                ov(i,j,k,1) = av(i,j,k,1) + ccy;
                ov(i,j,k,2) = av(i,j,k,2) + ccz;
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::RestrictLevel ( int a_lev )
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::RestrictLevel()");

    // nodal full weighting of the fine residual(s); fine ghosts are filled
    // (periodic/internal) or zero (physical: homogeneous Dirichlet). The
    // same index-space transfer serves the staggered components.
    auto restrict_fw = [&] (const MF& a_fine, MF& a_coarse)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a_coarse, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> const& fr = a_fine.const_array(mfi);
            Array4<Real> const& cb = a_coarse.array(mfi);
            amrex::ParallelFor(bx, 3,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
#if defined(WARPX_DIM_3D)
                    Real acc = 0.0_rt;
                    for (int oz = -1; oz <= 1; oz++) {
                    for (int oy = -1; oy <= 1; oy++) {
                    for (int ox = -1; ox <= 1; ox++) {
                        const Real w = (ox == 0 ? 0.5_rt : 0.25_rt)
                                     * (oy == 0 ? 0.5_rt : 0.25_rt)
                                     * (oz == 0 ? 0.5_rt : 0.25_rt);
                        acc += w * fr(2*i+ox, 2*j+oy, 2*k+oz, n);
                    }}}
                    cb(i,j,k,n) = acc;
#elif defined(WARPX_DIM_XZ)
                    Real acc = 0.0_rt;
                    for (int oz = -1; oz <= 1; oz++) {
                    for (int ox = -1; ox <= 1; ox++) {
                        const Real w = (ox == 0 ? 0.5_rt : 0.25_rt)
                                     * (oz == 0 ? 0.5_rt : 0.25_rt);
                        acc += w * fr(2*i+ox, 2*j+oz, k, n);
                    }}
                    cb(i,j,k,n) = acc;
#else // WARPX_DIM_1D_Z
                    cb(i,j,k,n) = 0.5_rt*fr(2*i,j,k,n)
                                + 0.25_rt*(fr(2*i-1,j,k,n) + fr(2*i+1,j,k,n));
#endif
                });
        }
    };

    restrict_fw(*m_rW[a_lev], *m_bW[a_lev+1]);
    if (m_pair) {
        restrict_fw(*m_rF[a_lev], *m_bF[a_lev+1]);
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ProlongLevel ( int a_lev )
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::ProlongLevel()");

    // The nodal (tri)linear prolongation reads coarse nodes ic and ic+1 <=
    // N/2 only -- never a coarse ghost -- so the exchange is needed only on
    // the cell-indexed staggered container (fine cell N-1 reads coarse
    // cell N/2, the periodic image of cell 0). The closure still applies
    // the wall-plane rows of the coarse correction. Level-0 index space
    // starts at 0, so i/2 floors correctly.
    auto prolong_lin = [&] (const MF& a_coarse, MF& a_fine)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> const& xc = a_coarse.const_array(mfi);
            Array4<Real> const& xf = a_fine.array(mfi);
            amrex::ParallelFor(bx, 3,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
#if defined(WARPX_DIM_3D)
                    const int ic = i/2; const int fx = i - 2*ic;
                    const int jc = j/2; const int fy = j - 2*jc;
                    const int kc = k/2; const int fz = k - 2*kc;
                    Real acc = 0.0_rt;
                    for (int c = 0; c <= fz; c++) {
                    for (int b = 0; b <= fy; b++) {
                    for (int a = 0; a <= fx; a++) {
                        acc += xc(ic+a, jc+b, kc+c, n);
                    }}}
                    xf(i,j,k,n) += acc * (fx ? 0.5_rt : 1.0_rt)
                                       * (fy ? 0.5_rt : 1.0_rt)
                                       * (fz ? 0.5_rt : 1.0_rt);
#elif defined(WARPX_DIM_XZ)
                    const int ic = i/2; const int fx = i - 2*ic;
                    const int jc = j/2; const int fy = j - 2*jc;
                    Real acc = 0.0_rt;
                    for (int b = 0; b <= fy; b++) {
                    for (int a = 0; a <= fx; a++) {
                        acc += xc(ic+a, jc+b, k, n);
                    }}
                    xf(i,j,k,n) += acc * (fx ? 0.5_rt : 1.0_rt)
                                       * (fy ? 0.5_rt : 1.0_rt);
#else // WARPX_DIM_1D_Z
                    const int ic = i/2; const int fx = i - 2*ic;
                    Real acc = 0.0_rt;
                    for (int a = 0; a <= fx; a++) { acc += xc(ic+a, j, k, n); }
                    xf(i,j,k,n) += acc * (fx ? 0.5_rt : 1.0_rt);
#endif
                });
        }
    };

    if (m_staggered) { FillGhosts(*m_xW[a_lev+1], a_lev+1); }
    ApplyDomainClosure(*m_xW[a_lev+1], a_lev+1, ClosureClass::JField);
    prolong_lin(*m_xW[a_lev+1], *m_xW[a_lev]);

    if (m_pair) {
        if (m_staggered) { FillGhosts(*m_xF[a_lev+1], a_lev+1); }
        ApplyDomainClosure(*m_xF[a_lev+1], a_lev+1, ClosureClass::EField);
        prolong_lin(*m_xF[a_lev+1], *m_xF[a_lev]);
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::Vcycle ( int a_lev, bool a_x_is_zero )
{
    if (a_lev == m_nlevels - 1) {
        const int nb = (m_nlevels == 1) ? m_nsweeps : kNuBottom;
        SmoothLevel(a_lev, nb, a_x_is_zero);
        return;
    }

    SmoothLevel(a_lev, kNuPre, a_x_is_zero);

    ResidualLevel(a_lev);
    FillGhosts(*m_rW[a_lev], a_lev);
    ApplyDomainClosure(*m_rW[a_lev], a_lev, ClosureClass::JField);
    if (m_pair) {
        FillGhosts(*m_rF[a_lev], a_lev);
        ApplyDomainClosure(*m_rF[a_lev], a_lev, ClosureClass::EField);
    }
    RestrictLevel(a_lev);
    if (!(m_zero_free && m_eb_flag == nullptr)) {
        m_xW[a_lev+1]->setVal(0.0);
        if (m_pair) { m_xF[a_lev+1]->setVal(0.0); }
    }
    Vcycle(a_lev + 1, true);
    ProlongLevel(a_lev);

    SmoothLevel(a_lev, kNuPost, false);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyMToW ( const MF& a_blk, int a_scomp )
{
    // block input -> m_bW[0] (+ m_bF[0]) in one kernel, then one V-cycle
    SplitBlock(a_blk, a_scomp);
    if (!(m_zero_free && m_eb_flag == nullptr)) {
        m_xW[0]->setVal(0.0);
        if (m_pair) { m_xF[0]->setVal(0.0); }
    }
    Vcycle(0, true);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::SplitBlock ( const MF& a_blk, int a_scomp )
{
    using namespace amrex;
    const bool two = (m_nkb == 6);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_bW[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real const> const& blk = a_blk.const_array(mfi, a_scomp);
        Array4<Real> const& bw = m_bW[0]->array(mfi);
        Array4<Real> const bf = two ? m_bF[0]->array(mfi) : Array4<Real>{};
        amrex::ParallelFor(bx, 3,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int c)
            {
                bw(i,j,k,c) = blk(i,j,k,c);
                if (two) { bf(i,j,k,c) = blk(i,j,k,3+c); }
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::MergeBlock ( MF& a_blk ) const
{
    using namespace amrex;
    const bool two = (m_nkb == 6);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_blk, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& blk = a_blk.array(mfi);
        Array4<Real const> const& xw = m_xW[0]->const_array(mfi);
        Array4<Real const> const xf = two ? m_xF[0]->const_array(mfi)
                                          : Array4<Real const>{};
        amrex::ParallelFor(bx, 3,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int c)
            {
                blk(i,j,k,c) = xw(i,j,k,c);
                if (two) { blk(i,j,k,3+c) = xf(i,j,k,c); }
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::BlockAxpy ( MF& a_dst, int a_dcomp, const MF* a_src, int a_nv,
                                     const amrex::GpuArray<amrex::Real, kMaxKry>& a_c,
                                     amrex::Real a_scale ) const
{
    using namespace amrex;
    const int nv = a_nv;
    const int nkb = m_nkb;
    const Real scale = a_scale;
    const bool have_src = (a_src != nullptr);
    const auto c = a_c;
    // Aliasing contract: the source is read at comps 0..nkb-1 and the
    // destination written at a_dcomp..a_dcomp+nkb-1, so in-place use is safe
    // only when both ranges coincide (a_dcomp == 0) or are disjoint; writing
    // into the Krylov block must not overlap the basis vectors being read.
    AMREX_ALWAYS_ASSERT(a_src != &a_dst || a_dcomp == 0 || a_dcomp >= nkb);
    AMREX_ALWAYS_ASSERT(&a_dst != m_kryB.get() || a_dcomp >= nv * nkb);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& d = a_dst.array(mfi, a_dcomp);
        Array4<Real const> const& V = m_kryB->const_array(mfi);
        Array4<Real const> const s = have_src ? a_src->const_array(mfi)
                                              : Array4<Real const>{};
        amrex::ParallelFor(bx, nkb,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int cc)
            {
                // in-place use (dst == src, a_dcomp == 0) is safe: each
                // thread then reads and writes the same element
                Real t = have_src ? s(i,j,k,cc) : 0.0_rt;
                for (int n = 0; n < nv; ++n) {
                    t -= c[n] * V(i,j,k, n*nkb + cc);
                }
                d(i,j,k,cc) = scale * t;
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::BatchDot ( const MF& a_w, int a_nv, bool a_with_ww, RT* a_out )
{
    using namespace amrex;
    using warpx::fields::FieldType;
    using ablastr::fields::Direction;

    const int nkb = m_nkb;
    const int nslots = a_nv + (a_with_ww ? 1 : 0);
    if (nslots == 0) { return; }
    // collocated grid: the nodal owner mask of Efield_fp (duplicated
    // periodic/box-shared nodes counted once); staggered containers are
    // cell-indexed and need no mask
    const bool use_mask = !m_staggered;
    const amrex::iMultiFab* mask = nullptr;
    if (use_mask) {
        mask = WarpX::GetInstance().getFieldDotMaskPointer(
            FieldType::Efield_fp, 0, Direction{0});
    }
    std::vector<amrex::Real> res(nslots, 0.0);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    if (Gpu::inLaunchRegion()) {
        constexpr int NV = 8;   // basis vectors per register chunk
        const int nbt = m_dot_nblocks_total;
        amrex::Real* part = m_dot_dev.data();
        for (MFIter mfi(a_w, false); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const int ib = mfi.LocalIndex();
            const int nblocks = m_dot_nblocks[ib];
            const int boff = m_dot_box_offset[ib];
            Array4<Real const> const w = a_w.const_array(mfi);
            Array4<Real const> const V = m_kryB->const_array(mfi);
            Array4<int const> const mk = use_mask ? mask->const_array(mfi)
                                                  : Array4<int const>{};
            const auto lo = amrex::lbound(bx);
            const auto len = amrex::length(bx);
            const Long npts = bx.numPts();
            const int nv = a_nv;
            const bool with_ww = a_with_ww;
            amrex::launch<kDotThreads>(nblocks, Gpu::gpuStream(),
                [=] AMREX_GPU_DEVICE () noexcept
                {
                    const int tid = static_cast<int>(threadIdx.x);
                    const Long stride = Long(blockDim.x) * Long(gridDim.x);
                    const Long start = Long(blockDim.x) * Long(blockIdx.x) + tid;
                    const Long lxy = Long(len.x) * Long(len.y);
                    // chunks of NV basis vectors; the <w,w> slot rides
                    // with the first chunk. Every thread executes every
                    // block reduction (uniform control flow).
                    const int nchunks = amrex::max(1, (nv + NV - 1) / NV);
                    for (int ch = 0; ch < nchunks; ++ch) {
                        const int n0 = ch * NV;
                        const bool do_ww = with_ww && (ch == 0);
                        Real s[NV];
                        for (int n = 0; n < NV; ++n) { s[n] = 0.0_rt; }
                        Real sww = 0.0_rt;
                        for (Long ic = start; ic < npts; ic += stride) {
                            const int k = static_cast<int>(ic / lxy);
                            const Long r = ic - Long(k) * lxy;
                            const int j = static_cast<int>(r / len.x);
                            const int i = static_cast<int>(r - Long(j) * len.x);
                            const int ii = i + lo.x;
                            const int jj = j + lo.y;
                            const int kk = k + lo.z;
                            if (use_mask && mk(ii,jj,kk) == 0) { continue; }
                            for (int c = 0; c < nkb; ++c) {
                                const Real wc = w(ii,jj,kk,c);
                                if (do_ww) { sww += wc * wc; }
                                for (int n = 0; n < NV; ++n) {
                                    if (n0 + n < nv) {
                                        s[n] += wc * V(ii,jj,kk,(n0 + n)*nkb + c);
                                    }
                                }
                            }
                        }
                        for (int n = 0; n < NV; ++n) {
                            const Real rs = Gpu::blockReduceSum<kDotThreads>(s[n]);
                            if (tid == 0 && n0 + n < nv) {
                                part[std::size_t(n0 + n) * nbt + boff + blockIdx.x] = rs;
                            }
                        }
                        if (do_ww) {
                            const Real rs = Gpu::blockReduceSum<kDotThreads>(sww);
                            if (tid == 0) {
                                part[std::size_t(nv) * nbt + boff + blockIdx.x] = rs;
                            }
                        }
                    }
                });
        }
        // the one host synchronization of the pass; fixed-order host sum
        const std::size_t ncopy = static_cast<std::size_t>(nslots) * nbt;
        Gpu::dtoh_memcpy(m_dot_host.data(), m_dot_dev.data(), ncopy * sizeof(amrex::Real));
        for (int n = 0; n < nslots; ++n) {
            amrex::Real ssum = 0.0_rt;
            for (int b = 0; b < nbt; ++b) {
                ssum += m_dot_host[static_cast<std::size_t>(n) * nbt + b];
            }
            res[n] = ssum;
        }
    } else
#endif
    {
        for (int n = 0; n < nslots; ++n) {
            const bool ww = (n == a_nv);
            const int vc = ww ? 0 : n * nkb;
            const MF& b = ww ? a_w : *m_kryB;
            res[n] = use_mask
                ? amrex::MultiFab::Dot(*mask, a_w, 0, b, vc, nkb, 0, true)
                : amrex::MultiFab::Dot(a_w, 0, b, vc, nkb, 0, true);
        }
    }
    if (amrex::ParallelContext::NProcsSub() > 1) {
        amrex::ParallelAllReduce::Sum(res.data(), nslots,
                                      amrex::ParallelContext::CommunicatorSub());
    }
    for (int n = 0; n < nslots; ++n) { a_out[n] = static_cast<RT>(res[n]); }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::RotateValid ( MF& a, int a_grow ) const
{
    using namespace amrex;
    // a <- (I+S)^-1 a pointwise (deposited-mass-matrix ion response,
    // see class docs); valid region grown by a_grow
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(a_grow);
        Array4<const Real> const& coef = m_coef[0]->const_array(mfi);
        Array4<Real> const& v = a.array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real vx = v(i,j,k,0);
                Real vy = v(i,j,k,1);
                Real vz = v(i,j,k,2);
                HybridPICPCInvIpS(coef(i,j,k,0), coef(i,j,k,1), coef(i,j,k,2),
                                  coef(i,j,k,6), coef(i,j,k,7), coef(i,j,k,8),
                                  vx, vy, vz);
                v(i,j,k,0) = vx;
                v(i,j,k,1) = vy;
                v(i,j,k,2) = vz;
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::RotateBInPlace ()
{
    RotateValid(*m_bE, 0);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::CopyCompIn ( const MF& a_src, MF& a_dst, int a_comp,
                                      bool a_add ) const
{
    using namespace amrex;
    // valid region of the container only; ghosts (including the top nodal
    // planes, which are ghost indices of the cell-indexed container) are
    // re-derived by FillBoundary and the closures
    const int nc = a_comp;
    const bool add = a_add;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<const Real> const& sv = a_src.const_array(mfi);
        Array4<Real> const& dv = a_dst.array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                dv(i,j,k,nc) = (add ? dv(i,j,k,nc) : 0.0_rt) + sv(i,j,k);
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::WriteXOut ( const MF& a_bE, const MF& a_t, T& a_x ) const
{
    using namespace amrex;
    for (int n = 0; n < 3; n++) {
        MF& x = *a_x.getArrayVec()[0][n];
        const int nc = n;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(x, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> const& bv = a_bE.const_array(mfi);
            Array4<const Real> const& tv = a_t.const_array(mfi);
            Array4<Real> const& xv = x.array(mfi);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    xv(i,j,k) = bv(i,j,k,nc) - tv(i,j,k,nc);
                });
        }
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::FillGhosts ( MF& a, int a_lev ) const
{
    using namespace amrex;

    const amrex::Geometry& geom = m_geom[a_lev];
    if (!m_fastfill_ok[a_lev]) {
        a.FillBoundary(geom.periodicity());
        return;
    }

    // Single domain-spanning box: every ghost cell that lies inside the
    // periodic domain image is the wrap of a valid cell. Index arithmetic
    // per dimension d with the CELL domain length N_d:
    //   in range [lo, hi]           -> unchanged (hi = lo+N nodal, lo+N-1 cell)
    //   out of range, periodic      -> lo + ((idx-lo) mod N)   in [lo, lo+N)
    //   out of range, not periodic  -> not written (closure territory)
    // A duplicated nodal plane index lo+N is mapped to lo so every ghost
    // derives from the primary plane. The valid region is never written.
    const int ncomp = a.nComp();
    GpuArray<int,3> per{0, 0, 0};
    GpuArray<int,3> len{1, 1, 1};
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        per[d] = geom.isPeriodic(d) ? 1 : 0;
        len[d] = geom.Domain().length(d);
    }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box vbx = mfi.validbox();
        // the grown tile (== the whole fab box without tiling, i.e. on GPU)
        const Box gbx = mfi.growntilebox();
        GpuArray<int,3> lo{0, 0, 0};
        GpuArray<int,3> hi{0, 0, 0};
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            lo[d] = vbx.smallEnd(d);
            hi[d] = vbx.bigEnd(d);
        }
        Array4<Real> const& v = a.array(mfi);
        amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int s[3] = {i, j, k};
                bool ghost = false;
                for (int d = 0; d < AMREX_SPACEDIM; d++) {
                    if (s[d] < lo[d] || s[d] > hi[d]) {
                        if (!per[d]) { return; }
                        ghost = true;
                        int r = (s[d] - lo[d]) % len[d];
                        if (r < 0) { r += len[d]; }
                        s[d] = lo[d] + r;
                    } else if (per[d] && s[d] == lo[d] + len[d]) {
                        s[d] = lo[d];   // duplicated nodal plane -> primary
                    }
                }
                if (!ghost) { return; }
                for (int n = 0; n < ncomp; n++) {
                    v(i,j,k,n) = v(s[0],s[1],s[2],n);
                }
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyDomainClosure ( MF& a, int a_lev,
                                              ClosureClass a_cls,
                                              int a_scomp, int a_ncomp ) const
{
    using namespace amrex;

    // component range: the closure logic indexes components RELATIVE to
    // a_scomp (so a 6-comp Krylov block gets its W-half and E-half closed
    // with their own classes by two calls), the storage by a_scomp + n
    const int ncl = (a_ncomp < 0) ? a.nComp() : a_ncomp;
    const int sc = a_scomp;

    if (m_has_pec_wall) {
        const Box dom = amrex::convert(m_geom[a_lev].Domain(),
                                       a.boxArray().ixType());
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            for (int side = 0; side < 2; side++) {
                if (!m_wall_pec[d][side]) { continue; }
                // wall NODE index: the staggered containers are
                // cell-indexed, so the hi-side wall plane sits one past
                // the cell bigEnd
                const int w = (side == 0)
                    ? dom.smallEnd(d)
                    : dom.bigEnd(d) + (m_staggered ? 1 : 0);
                const int nc = m_dim2comp[d];
                const bool lo_side = (side == 0);
                const int cls = static_cast<int>(a_cls);
                const int dd = d;
                const bool stag = m_staggered;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.growntilebox();
                    // component window: relative index n, storage sc + n
                    Array4<Real> const& v = a.array(mfi, sc);
                    amrex::ParallelFor(bx, ncl,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            const int id = (dd == 0) ? i : ((dd == 1) ? j : k);
                            const bool outside =
                                lo_side ? (id < w) : (id > w);
                            // Coef: comps 0-2 carry bhat (vector), the rest
                            // are scalars
                            const bool scalar =
                                (cls == 4/*Coef*/) && (n >= 3);
                            const bool is_norm = !scalar && (n == nc);
                            if (stag) {
                                // Yee: components NODAL along dd (wall-
                                // tangential E/J, wall-normal B) carry a
                                // wall-plane DOF that is zero at PEC (dead
                                // tangential-E columns, zero normal B, and
                                // zero mirror-image tangential J) and
                                // mirror ODD about the wall node;
                                // components CELL-centered along dd have no
                                // wall DOF and mirror EVEN about the wall
                                // plane. The node-built coefficients keep a
                                // smooth even cell-mirror (their ghost ring
                                // only shapes U there).
                                bool nodal_d;
                                Real sgn;
                                if (cls == 0/*EField*/ || cls == 2/*JField*/) {
                                    nodal_d = !is_norm;
                                    sgn = is_norm ? 1.0_rt : -1.0_rt;
                                } else if (cls == 4/*Coef*/) {
                                    nodal_d = false;
                                    sgn = 1.0_rt;
                                } else { // BField, BCoef
                                    nodal_d = is_norm;
                                    sgn = is_norm ? -1.0_rt : 1.0_rt;
                                }
                                int im;
                                if (nodal_d) {
                                    const bool beyond =
                                        lo_side ? (id < w) : (id > w);
                                    if (!beyond) {
                                        if (id == w) { v(i,j,k,n) = 0.0_rt; }
                                        return;
                                    }
                                    im = 2*w - id;
                                } else {
                                    const bool beyond =
                                        lo_side ? (id < w) : (id >= w);
                                    if (!beyond) { return; }
                                    im = 2*w - 1 - id;
                                }
                                const int mi = (dd == 0) ? im : i;
                                const int mj = (dd == 1) ? im : j;
                                const int mk = (dd == 2) ? im : k;
                                v(i,j,k,n) = sgn * v(mi,mj,mk,n);
                                return;
                            }
                            if (!outside) {
                                if (id != w) { return; }
                                // wall plane: EField zeroes tangential
                                // (identity rows / dead columns), BField
                                // zeroes the normal component
                                if (cls == 0/*EField*/ && !scalar
                                    && !is_norm) { v(i,j,k,n) = 0.0_rt; }
                                if (cls == 1/*BField*/ && is_norm) {
                                    v(i,j,k,n) = 0.0_rt;
                                }
                                return;
                            }
                            // ghost planes: signed mirror about the wall node
                            const int im = 2*w - id;
                            Real sign = 1.0_rt;
                            if (!scalar) {
                                if (cls == 0/*EField*/ || cls == 2/*JField*/) {
                                    sign = is_norm ? 1.0_rt : -1.0_rt;
                                } else { // BField, BCoef, Coef bhat comps
                                    sign = is_norm ? -1.0_rt : 1.0_rt;
                                }
                            }
                            const int mi = (dd == 0) ? im : i;
                            const int mj = (dd == 1) ? im : j;
                            const int mk = (dd == 2) ? im : k;
                            v(i,j,k,n) = sign * v(mi,mj,mk,n);
                        });
                }
            }
        }
    }

#if defined(WARPX_DIM_RZ)
    // Symmetry axis (r = 0): signed-mirror guards about the axis node --
    // r/theta components odd, z even, the same parities for E-, B- and
    // J-like objects at m = 0 (matching ApplyFieldBoundaryOnAxis). The Ohm
    // solve zeroes E_theta ON the axis, so EField-class objects get that
    // row zeroed here as well (restored as an identity row on the output).
    if (m_has_axis) {
        const Box adom = amrex::convert(m_geom[a_lev].Domain(),
                                        a.boxArray().ixType());
        const int w = adom.smallEnd(0);
        const int cls = static_cast<int>(a_cls);
        const bool stag = m_staggered;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox();
            Array4<Real> const& v = a.array(mfi, sc);
            amrex::ParallelFor(bx, ncl,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    if (i > w) { return; }
                    const bool scalar = (cls == 4/*Coef*/) && (n >= 3);
                    // Yee: components CELL-centered in r (the E/J r
                    // component, the B theta/z components; never the
                    // node-built coefficients) have no axis DOF and mirror
                    // with the residual's ishift = 1
                    bool cell_r = false;
                    if (stag && cls != 4/*Coef*/) {
                        cell_r = (cls == 0 || cls == 2) ? (n == 0)
                                                        : (n != 0);
                    }
                    if (i == w) {
                        if (cls == 0/*EField*/ && n == 1) {
                            v(i,j,k,n) = 0.0_rt;
                        }
                        return;
                    }
                    const int im = cell_r ? (2*w - 1 - i) : (2*w - i);
                    const Real sign =
                        (scalar || n == 2) ? 1.0_rt : -1.0_rt;
                    v(i,j,k,n) = sign * v(im,j,k,n);
                });
        }
    }
#endif

    // EB (staircase): zero covered/frozen nodes so no stencil emits or
    // receives contributions there (the residual's masked Ohm solves skip
    // those locations). Fine level only; masked up to the flag's guards.
    // Exception: the pair's local W-hat scale (Coef comp 14) is a strictly
    // positive SCALING, divided by in the joint block -- zeroing it at
    // covered nodes makes 1/sW = inf there (NaN poisoning that the
    // level-0-only closure cannot clean off the coarse levels).
    if (a_lev == 0 && m_eb_flag != nullptr) {
        amrex::IntVect ngm = a.nGrowVect();
        ngm.min(m_eb_flag->nGrowVect());
        const bool is_coef = (a_cls == ClosureClass::Coef);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(ngm);
            Array4<Real> const& v = a.array(mfi, sc);
            Array4<const int> const& fl = m_eb_flag->const_array(mfi);
            amrex::ParallelFor(bx, ncl,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    if (is_coef && n == 14) { return; }
                    if (fl(i,j,k) == 0) { v(i,j,k,n) = 0.0_rt; }
                });
        }
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::OverrideIdentityRows ( T& a_out, const T& a_in ) const
{
    using namespace amrex;

    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        for (int side = 0; side < 2; side++) {
            // PEC faces: exact identity tangential rows (both grids).
            // Staggered hi-side zero-extension faces: the nodal wall
            // plane is an unmodeled ghost plane of the cell-indexed
            // container -- the preconditioner would return zero there
            // (a singular preconditioner); make those rows identity.
            const bool open_stag = m_staggered && (side == 1)
                && m_wall_open[d][side];
            if (!m_wall_pec[d][side] && !open_stag) { continue; }
            const int nc = m_dim2comp[d];
            for (int n = 0; n < 3; n++) {
                if (n == nc) { continue; }  // tangential rows only
                MF& o = *a_out.getArrayVec()[0][n];
                const MF& s = *a_in.getArrayVec()[0][n];
                // wall plane in this component's own staggering (the
                // tangential components are nodal along d on both grids)
                const Box dom = amrex::convert(m_geom[0].Domain(),
                                               o.boxArray().ixType());
                const int w = (side == 0) ? dom.smallEnd(d) : dom.bigEnd(d);
                Box plane = dom;
                plane.setSmall(d, w);
                plane.setBig(d, w);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(o, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box bx = mfi.tilebox() & plane;
                    if (bx.isEmpty()) { continue; }
                    Array4<Real> const& ov = o.array(mfi);
                    Array4<const Real> const& sv = s.const_array(mfi);
                    amrex::ParallelFor(bx,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        { ov(i,j,k) = sv(i,j,k); });
                }
            }
        }
    }

#if defined(WARPX_DIM_RZ)
    // E_theta on the axis is fixed to zero by the Ohm solve: identity row
    if (m_has_axis) {
        const Box adom = amrex::convert(m_geom[0].Domain(),
            a_out.getArrayVec()[0][1]->boxArray().ixType());
        Box plane = adom;
        plane.setBig(0, adom.smallEnd(0));
        MF& o = *a_out.getArrayVec()[0][1];
        const MF& si = *a_in.getArrayVec()[0][1];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(o, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box bx = mfi.tilebox() & plane;
            if (bx.isEmpty()) { continue; }
            Array4<Real> const& ov = o.array(mfi);
            Array4<const Real> const& sv = si.const_array(mfi);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                { ov(i,j,k) = sv(i,j,k); });
        }
    }
#endif

    if (m_eb_flag != nullptr) {
        for (int n = 0; n < 3; n++) {
            MF& o = *a_out.getArrayVec()[0][n];
            const MF& s = *a_in.getArrayVec()[0][n];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(o, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                Array4<Real> const& ov = o.array(mfi);
                Array4<const Real> const& sv = s.const_array(mfi);
                Array4<const int> const& fl = m_eb_flag->const_array(mfi);
                amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (fl(i,j,k) == 0) { ov(i,j,k) = sv(i,j,k); }
                    });
            }
        }
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::InnerGMRES ()
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::InnerGMRES()");

    const int m = m_inner_max;
    const int nkb = m_nkb;
    // pair: the s_W-weighted block norm is E-space accurate, so inner_rtol
    // applies directly. W-system: the E-space residual of the reconstructed
    // solution is r_E = -U r_W, so the W-space tolerance is divided by the
    // grid-scale amplification of U (up to beta_max).
    const RT rtol = m_pair
        ? m_inner_rtol
        : m_inner_rtol / amrex::max(RT(1.0), m_beta_max);

    // v_0 = b/||b|| (the RHS block was assembled by Apply() into m_bKB)
    amrex::MultiFab::Copy(*m_kryB, *m_bKB, 0, 0, nkb, 0);
    m_xW[0]->setVal(0.0);
    if (nkb == 6) { m_xF[0]->setVal(0.0); }
    std::vector<RT> d1(m + 2, RT(0.0));
    std::vector<RT> d2(m + 2, RT(0.0));
    BatchDot(*m_kryB, 0, true, d1.data());
    const RT beta = std::sqrt(d1[0]);
    if (beta <= RT(0.0)) { return; }
    m_kryB->mult(RT(1.0)/beta, 0, nkb, 0);

    std::vector<std::vector<RT>> H(m + 1, std::vector<RT>(m, RT(0.0)));
    std::vector<RT> cs(m, RT(0.0));
    std::vector<RT> sn(m, RT(0.0));
    std::vector<RT> g(m + 1, RT(0.0));
    g[0] = beta;
    amrex::GpuArray<amrex::Real, kMaxKry> cb{};

    int its = 0;
    for (int j = 0; j < m; j++) {
        its = j + 1;

        // z = M v_j (one V-cycle; result merged into the ghost-filled block
        // m_zB), w = A z
        ApplyMToW(*m_kryB, j*nkb);
        MergeBlock(*m_zB);
        FillGhosts(*m_zB, 0);
        ApplyDomainClosure(*m_zB, 0, ClosureClass::JField, 0, 3);
        if (nkb == 6) {
            ApplyDomainClosure(*m_zB, 0, ClosureClass::EField, 3, 3);
        }
        ApplyReducedOp(*m_zB, *m_wB);

        // Classical Gram-Schmidt with reorthogonalization (CGS2):
        //   pass 1: d1 = V^T w,        w1 = w  - V d1
        //   pass 2: d2 = V^T w1,       w2 = w1 - V d2,  H_j = d1 + d2,
        //   h_{j+1,j} = ||w2|| = sqrt(||w1||^2 - ||d2||^2)  (Pythagoras; an
        //   explicit norm is taken when the second projection is not tiny)
        // Two batched reductions (two host syncs) per iteration instead of
        // the j+2 sequential ones of modified Gram-Schmidt; identical in
        // exact arithmetic.
        BatchDot(*m_wB, j + 1, false, d1.data());
        for (int i = 0; i <= j; i++) { cb[i] = static_cast<amrex::Real>(d1[i]); }
        BlockAxpy(*m_wB, 0, m_wB.get(), j + 1, cb, 1.0_rt);
        BatchDot(*m_wB, j + 1, true, d2.data());
        RT d2sq = RT(0.0);
        for (int i = 0; i <= j; i++) {
            H[i][j] = d1[i] + d2[i];
            d2sq += d2[i]*d2[i];
        }
        RT h2 = d2[j+1] - d2sq;
        const bool explicit_norm = (d2sq > RT(1.0e-8) * d2[j+1]) || (h2 <= RT(0.0));
        if (explicit_norm) {
            // near-breakdown regime: form w2 and measure it directly
            for (int i = 0; i <= j; i++) { cb[i] = static_cast<amrex::Real>(d2[i]); }
            BlockAxpy(*m_wB, 0, m_wB.get(), j + 1, cb, 1.0_rt);
            RT ww[1] = {RT(0.0)};
            BatchDot(*m_wB, 0, true, ww);
            h2 = ww[0];
        }
        const RT hj1 = std::sqrt(amrex::max(h2, RT(0.0)));

        for (int i = 0; i < j; i++) {
            const RT t = cs[i]*H[i][j] + sn[i]*H[i+1][j];
            H[i+1][j] = -sn[i]*H[i][j] + cs[i]*H[i+1][j];
            H[i][j] = t;
        }
        const RT denom = std::sqrt(H[j][j]*H[j][j] + hj1*hj1);
        cs[j] = H[j][j] / denom;
        sn[j] = hj1 / denom;
        H[j][j] = denom;
        g[j+1] = -sn[j]*g[j];
        g[j] = cs[j]*g[j];

        const bool converged = (std::abs(g[j+1]) <= rtol * beta);
        const bool breakdown = (hj1 <= RT(1.0e-14) * beta);
        if (converged || breakdown || j + 1 == m) { break; }

        // v_{j+1} = w2 / h_{j+1,j}, fused with the second projection when
        // w2 was not materialized above
        if (explicit_norm) {
            BlockAxpy(*m_kryB, (j+1)*nkb, m_wB.get(), 0, cb, 1.0_rt/static_cast<amrex::Real>(hj1));
        } else {
            for (int i = 0; i <= j; i++) { cb[i] = static_cast<amrex::Real>(d2[i]); }
            BlockAxpy(*m_kryB, (j+1)*nkb, m_wB.get(), j + 1, cb,
                      1.0_rt/static_cast<amrex::Real>(hj1));
        }
    }

    std::vector<RT> y(its, RT(0.0));
    for (int i = its - 1; i >= 0; i--) {
        RT s = g[i];
        for (int k = i + 1; k < its; k++) { s -= H[i][k]*y[k]; }
        y[i] = s / H[i][i];
    }

    // x = M (V y): the V-cycle is linear, so apply it once to the summed
    // basis (classical right-preconditioned GMRES update); one fused
    // linear-combination kernel forms V y (BlockAxpy with no source term
    // computes -sum c_i v_i, hence the sign)
    for (int i = 0; i < its; i++) { cb[i] = -static_cast<amrex::Real>(y[i]); }
    BlockAxpy(*m_zB, 0, nullptr, its, cb, 1.0_rt);
    ApplyMToW(*m_zB, 0);   // final W (and E) in m_xW[0]/m_xF[0]

    if (m_verbose) {
        amrex::Print() << "  pc_hybrid_pic inner GMRES: " << its << " its, rel res "
                       << std::abs(g[its])/beta << "\n";
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::Apply ( T& a_x, const T& a_b )
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::Apply()");

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(), "HybridPICPC::Apply() called on undefined object");

    // pe-as-Newton-unknown triangular leg (see m_pe_row_active): solve the
    // pe row with the pointwise diagonal D_p, then fold its Ohm gradient
    // coupling into the E-row RHS below (block upper-triangular approximate
    // Jacobian; the reverse pe <- E coupling is dropped). With the leg
    // inactive, or an E-only vector, the pe row falls through to the
    // identity passthrough at the end of this function.
    const bool pe_row = m_pe_row_active
        && (a_x.getScalarVecType() != warpx::fields::FieldType::None);
    if (pe_row) {
        amrex::MultiFab& xp = *a_x.getScalarVec()[0];
        const amrex::MultiFab& bp = *a_b.getScalarVec()[0];
        // the scheme scales the pe Newton row; the E leg needs physical dpe
        const amrex::Real inv_kappa =
            1.0_rt / m_ops->GetHybridPICModel()->m_pe_newton_scale;
        m_pe_row_scr->setVal(0.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(xp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box vbx = mfi.tilebox();
            auto const& x  = xp.array(mfi);
            auto const& b  = bp.const_array(mfi);
            auto const& cf = m_coef[0]->const_array(mfi);
            auto const& s  = m_pe_row_scr->array(mfi);
            amrex::ParallelFor(vbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    const Real v = b(i,j,k) / (1.0_rt + cf(i,j,k,10));
                    x(i,j,k) = v;
                    s(i,j,k) = v * inv_kappa;   // physical dpe for the E leg
                });
        }
        FillGhosts(*m_pe_row_scr, 0);
    }

    // b~: the PC input b (rotated by (I+S)^-1 with the ion response), with
    // the pe-row gradient folded in and the E-field closure applied
    for (int n = 0; n < 3; n++) {
        if (m_staggered) {
            CopyCompIn(*a_b.getArrayVec()[0][n], *m_bE, n, false);
        } else {
            amrex::MultiFab::Copy(*m_bE, *a_b.getArrayVec()[0][n], 0, n, 1, 0);
        }
    }
    if (pe_row) {
        // b_E <- b_E - grad(dpe)/rho_e: the dF_E/dpe Ohm coupling the
        // identity pe leg left entirely to GMRES. Domain ghosts of the dpe
        // scratch are zero-extended; the affected wall-tangential E rows are
        // identity-overridden on output anyway.
        const auto dxi = m_geom[0].InvCellSizeArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*m_bE, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box vbx = mfi.tilebox();
            auto const& be = m_bE->array(mfi);
            auto const& s  = m_pe_row_scr->const_array(mfi);
            auto const& cf = m_coef[0]->const_array(mfi);
            amrex::ParallelFor(vbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    const Real ir = 1.0_rt / cf(i,j,k,9);
#if defined(WARPX_DIM_1D_Z)
                    be(i,j,k,2) -= ir * (s(i+1,j,k) - s(i-1,j,k)) * 0.5_rt * dxi[0];
#elif defined(WARPX_DIM_XZ)
                    be(i,j,k,0) -= ir * (s(i+1,j,k) - s(i-1,j,k)) * 0.5_rt * dxi[0];
                    be(i,j,k,2) -= ir * (s(i,j+1,k) - s(i,j-1,k)) * 0.5_rt * dxi[1];
#elif defined(WARPX_DIM_3D)
                    be(i,j,k,0) -= ir * (s(i+1,j,k) - s(i-1,j,k)) * 0.5_rt * dxi[0];
                    be(i,j,k,1) -= ir * (s(i,j+1,k) - s(i,j-1,k)) * 0.5_rt * dxi[1];
                    be(i,j,k,2) -= ir * (s(i,j,k+1) - s(i,j,k-1)) * 0.5_rt * dxi[2];
#else
                    amrex::ignore_unused(i, j, k, ir, be, s, dxi);
#endif
                });
        }
    }
    if (m_have_mm) {
        // the W-system RHS is K (I+S)^-1 b: rotate b in place first (the
        // reconstruction below re-reads m_bE, which must stay the ROTATED
        // b so that x = (I+S)^-1 b - (I+S)^-1 U W)
        RotateBInPlace();
    }
    FillGhosts(*m_bE, 0);
    // PEC dead tangential-E columns (wall values zeroed; the identity rows
    // are restored on the output below) + parity ghosts; EB nodes zeroed
    ApplyDomainClosure(*m_bE, 0, ClosureClass::EField);
    m_xW[0]->setVal(0.0);
    if (m_pair) {
        // primal pair: E-row RHS (block comps 3..5) is b~ itself, W-row RHS
        // (comps 0..2) is zero
        m_bKB->setVal(0.0, 0, 3, 0);
        amrex::MultiFab::Copy(*m_bKB, *m_bE, 0, 3, 3, 0);
    } else {
        // bW = curl curl b~ (exact operators, one exchange): the curl-free
        // content of b never enters the W-solve and passes through the
        // reconstruction exactly
        ComputeCurl(*m_bE, *m_tmp[0], 0, 1, CurlKind::EtoB);
        ApplyDomainClosure(*m_tmp[0], 0, ClosureClass::BField);
        MF& bdst = m_nested ? *m_bKB : *m_bW[0];
        ComputeCurl(*m_tmp[0], bdst, 0, 0, CurlKind::BtoJ);
    }

    if (m_nested) {
        InnerGMRES();
    } else {
        SmoothLevel(0, m_nsweeps, true);
    }

    if (m_pair) {
        // primal pair: the E-half m_xF[0] IS the output; route it through
        // the common x = bE - tmp emitter below
        amrex::MultiFab::LinComb(*m_tmp[0], 1.0_rt, *m_bE, 0,
                                 -1.0_rt, *m_xF[0], 0, 0, 3, 0);
        if (m_staggered) {
            // WriteXOut reads the top nodal planes from the first ghost ring
            // of the cell container: fill it (periodic pair grids only)
            FillGhosts(*m_tmp[0], 0);
        }
    } else {
        // Reconstruct x = (I+S)^-1 (b - U xW) = bE_rot - (I+S)^-1 U xW
        // (m_bE already holds the rotated b when the ion-response leg is on)
        FillGhosts(*m_xW[0], 0);
        ApplyDomainClosure(*m_xW[0], 0, ClosureClass::JField);
        // staggered: evaluate one ghost ring so the top nodal planes of the
        // reconstruction can be read back by WriteXOut
        const int gout = m_staggered ? 1 : 0;
        ApplyUValid(*m_xW[0], *m_tmp[0], 0, gout);
        if (m_have_mm) { RotateValid(*m_tmp[0], gout); }
        // the residual BCs the Ohm output: wall-tangential/EB contributions off
        ApplyDomainClosure(*m_tmp[0], 0, ClosureClass::EField);
    }
    if (m_staggered) {
        WriteXOut(*m_bE, *m_tmp[0], a_x);
    } else {
        for (int n = 0; n < 3; n++) {
            amrex::MultiFab::LinComb(*a_x.getArrayVec()[0][n],
                                     1.0_rt, *m_bE, n,
                                     -1.0_rt, *m_tmp[0], n,
                                     0, 1, 0);
        }
    }
    // exact identity rows of the true Jacobian: x = b on wall-tangential E
    // at PEC faces and on EB-covered nodes
    OverrideIdentityRows(a_x, a_b);

    // pe-as-Newton-unknown fallback: with the triangular pe-row leg
    // inactive, pass the scalar (pe) row through as identity so the
    // right-preconditioned operator is nonsingular in those rows (the
    // precond entry zeroes a_x first, so skipping this would pin dpe = 0).
    // With the leg active, a_x's scalar already holds the D_p solve above.
    if (!pe_row && a_x.getScalarVecType() != warpx::fields::FieldType::None) {
        amrex::MultiFab::Copy(*a_x.getScalarVec()[0], *a_b.getScalarVec()[0],
                              0, 0, 1, 0);
    }
}

#endif
