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
 *  Because the dE coefficient is the identity (electron convection dropped),
 *  dE is eliminated EXACTLY in favor of W = curl curl dE -- the
 *  parabolization of Chacon 2025 Sec. 3 applied globally (the A-based row
 *  used there does not admit it):
 *
 *    (I + curlcurl U) W = curlcurl b,      dE = b - U W.
 *
 *  The reduction is essential for an E-based solver: the curl-free
 *  (electrostatic) content of b -- physical and large in hybrid PIC, unlike
 *  the gauge content of A -- passes through the reconstruction exactly
 *  instead of being 1/beta-suppressed by the defect-corrected block diagonal
 *  (beta = l2 c_w, the whistler CFL at the grid Nyquist; l2 = sum_d 2/dx_d^2).
 *
 *  The W-system is smoothed with damped block-Jacobi (Chacon 2025 Sec. 4.1): the
 *  defect correction builds the per-node diagonal with curlcurl ~ -del2 while
 *  all residuals use the exact central-difference operators, so the fixed
 *  point is exact. The diagonal P = alpha I + beta (. x bhat), alpha =
 *  1 + l2 (a_e + a_H l2), has the analytic inverse
 *
 *    P^-1 v = [alpha v - beta (v x bhat) + (beta^2/alpha)(v.bhat) bhat]
 *             / (alpha^2 + beta^2).
 *
 *  Modes (pc_hybrid_pic.mode):
 *
 *  "smoother" (default): nsweeps damped Jacobi sweeps on the fine grid.
 *    Convergence of the stiffest regions needs nsweeps ~ 1.5-3x the largest
 *    beta (grows linearly with dt); sweeps are field-only and ~10x cheaper
 *    than a Krylov matvec, so this is the robust default under the plain
 *    (non-flexible) amrex_gmres outer solver, where the preconditioner must
 *    be a fixed linear operator.
 *
 *  "nested": the Sec.-4 structure of Chacon 2025. Geometric-MG V(3,3) (factor-2
 *    coarsening while the boxes allow, same DistributionMapping; coefficient
 *    injection, l2 rebuilt per level; nodal full-weighting restriction,
 *    (tri)linear prolongation; 8 bottom sweeps) preconditions an inner
 *    right-preconditioned GMRES on the W-system. The Krylov space absorbs
 *    the per-level Nyquist rings of the central-difference curl that
 *    smoothing alone cannot. Convergence targets E-SPACE accuracy: from
 *    dE = b - U W, the E-residual of the reconstructed solution is exactly
 *    r_E = -U r_W, so a W-space tolerance overstates convergence by up to
 *    the amplification of U (~beta at the grid scale). The inner W-space
 *    tolerance is therefore INNER_RTOL / max(1, beta_max), with beta_max =
 *    l2 max(c_w) refrozen each Newton iteration, and up to INNER_MAX
 *    iterations (the reference 10-iteration/0.1 recipe is recovered when
 *    beta_max <= 1). The inner solve makes the preconditioner NONLINEAR: it
 *    requires a flexible outer Krylov method (newton.linear_solver =
 *    petsc_ksp with PETSC_OPTIONS="-ksp_type fgmres").
 *
 *  Performance layout (launch-latency optimizations; all multi-rank safe):
 *  the three W components live in single 3-component MultiFabs, so every
 *  ghost exchange is one FillBoundary call (one message round) instead of
 *  three; iterates carry 3 ghost layers and one exchange feeds the whole
 *  U -> curl -> curl composition per sweep, with U and the first curl
 *  evaluated into ghost regions. One sweep = 1 FillBoundary + 3 kernels.
 *
 *  Physical boundaries: PEC(_Insulator) faces are modeled exactly
 *  (identity tangential-E rows with dead columns and parity-mirrored
 *  delta-B/delta-J ghosts, matching the residual's boundary application;
 *  see ClosureClass); other non-periodic BCs fall back to a homogeneous
 *  zero-extension closure. With embedded boundaries, EB-covered nodes are
 *  masked off the same nodal eb_update_E flag the residual's Ohm solves
 *  consume and act as identity rows.
 *
 *  With implicit_evolve.use_mass_matrices_pc = 1 the deposited reduced
 *  mass-matrix diagonal Sigma = dJ/dE (the ion response, as consumed by the
 *  EM-PIC preconditioners) is folded through the Hall term as
 *  J_pc = (I + S) + U K, S v = (Sigma v) x B0 / rho_c, with the exact
 *  reduction preserved via the pointwise closed-form (I+S)^-1.
 *
 *  Geometries: collocated (nodal) Cartesian 1D/2D/3D, and RZ restricted to
 *  the m = 0 azimuthal mode (the residual's CylindricalYee up/down stencil
 *  compositions, on-axis branches and axis mirror parities are reproduced
 *  exactly; the r = 0 axis and PEC outer wall are part of the boundary
 *  closure).
 *  Default damping sigma = 0.55 (the 3D Nyquist stability bound of the
 *  damped Jacobi sweep is sigma < 2/3).
 *
 *  Input parameters (prefix pc_hybrid_pic.):
 *    mode     (string, default "smoother")  "smoother" | "nested"
 *    nsweeps  (int,    default 48)          sweeps in smoother mode
 *    sigma    (float,  default 0.55)        Jacobi damping, in (0, 2)
 *    mg_floor (int,    default 4)           coarsest MG box width (nested)
 *    verbose  (bool,   default false)
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

        /** Approximately solve J_pc x = b via the reduced W-system. */
        void Apply (T& a_x, const T& a_b) override;

        void printParameters () const override;

        [[nodiscard]] inline bool IsDefined () const override { return m_is_defined; }

        /* Implementation details; public because CUDA requires the enclosing
         * function of an extended device lambda to have public access. */

        /** tmp = U a on growntilebox(2); tmp2 = curl tmp on growntilebox(1).
         *  a must be ghost-filled (3 layers). After this, curl(tmp2) is
         *  available on the valid region without further exchanges. */
        void StageKU (const MF& a, int a_lev);

        /** Fused damped block-Jacobi sweeps on level a_lev: per sweep, one
         *  FillBoundary of xW then three kernels (StageKU + update). */
        void SmoothLevel (int a_lev, int a_nsweeps, bool a_x_is_zero);

        /** rW = bW - xW - curl curl (U xW) on the valid region of a_lev. */
        void ResidualLevel (int a_lev);

        /** out = a + curl curl (U a) on the fine level (a ghost-filled). */
        void ApplyReducedOp (const MF& a, MF& out);

        /** Full-weighting restriction of rW[a_lev] into bW[a_lev+1]. */
        void RestrictLevel (int a_lev);

        /** (Tri)linear prolongation of xW[a_lev+1] added into xW[a_lev]. */
        void ProlongLevel (int a_lev);

        /** Coefficient injection from level a_lev-1 to a_lev. */
        void RestrictCoef (int a_lev);

        /** out = curl a (valid region grown by a_grow; a must have ghosts). */
        /** Which discrete curl composition stage this is: the Faraday curl
         *  (E-like -> B-like, upward stencils in RZ) or the Ampere curl
         *  (B-like -> J-like, downward stencils with the on-axis branches
         *  in RZ). Cartesian uses the same central stencil for both. */
        enum class CurlKind { EtoB, BtoJ };

        void ComputeCurl (const MF& a, MF& out, int a_lev, int a_grow,
                          CurlKind a_kind) const;

        /** Final reconstruction helper: out = U a on the valid region. */
        void ApplyUValid (const MF& a, MF& out, int a_lev) const;

        /** Pointwise (I+S)^-1 on the valid region (ion-response leg). */
        void RotateValid (MF& a) const;
        void RotateBInPlace ();

        /** Physical-boundary closure classes for the vector objects in the
         *  PC pipeline, mirroring the residual's PEC handling exactly
         *  (see ApplyDomainClosure): at a PEC face with normal component n,
         *   EField : linearized E-like objects (the PC input b and the
         *            linearized Ohm outputs U(.)): tangential components
         *            zeroed ON the wall plane (ApplyEfieldBoundary is applied
         *            to the iterate AND to the Ohm output in the residual),
         *            ghosts mirror-filled tangential-odd / normal-even;
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
         *            signs, comps 3-8 even; mirror ghosts, no wall write). */
        enum class ClosureClass { EField, BField, JField, BCoef, Coef };

        /** Apply the physical-boundary closure of class a_cls to a on level
         *  a_lev (PEC/PEC_Insulator faces only; other non-periodic BCs keep
         *  the previous zero-extension), then zero EB-covered nodes on the
         *  fine level when an EB mask is present. */
        void ApplyDomainClosure (MF& a, int a_lev, ClosureClass a_cls) const;

        /** Enforce the exact identity rows of the true Jacobian on a_out:
         *  wall-tangential components at PEC faces get a_out = a_in, and
         *  EB-covered nodes get a_out = a_in (the PC acts as the identity
         *  there; the covered residual rows are structurally zero). */
        void OverrideIdentityRows (T& a_out, const T& a_in) const;

    private:

        /** Recursive V(3,3) cycle for the correction on level a_lev. */
        void Vcycle (int a_lev, bool a_x_is_zero);

        /** One preconditioner application for the inner GMRES (one V-cycle;
         *  nsweeps sweeps if the grid is not coarsenable); result in m_xW[0]. */
        void ApplyMToW (const MF& a_in);

        /** Masked (owned-node) dot product over the three components. */
        [[nodiscard]] RT DotW (const MF& a, const MF& b) const;

        /** Inner GMRES on the W-system (<= INNER_ITERS its or INNER_RTOL);
         *  result in m_xW[0]. */
        void InnerGMRES ();

        // Paper-faithful fixed constants (Sec. 4: V(3,3) smoothing; EMHD
        // block inverted by MG-preconditioned GMRES, <= 10 its or rtol 0.1)
        static constexpr int NU_PRE = 3;
        static constexpr int NU_POST = 3;
        static constexpr int NU_BOTTOM = 8;
        static constexpr int INNER_MAX = 60;
        static constexpr double INNER_RTOL = 0.1;
        RT m_sigma = 0.55;
        /** Ion-response leg active: the deposited reduced mass-matrix
         *  diagonal (raw dJ/dE) is available via
         *  implicit_evolve.use_mass_matrices_pc = 1 (set in Update()). */
        bool m_have_mm = false;
        /** Iterate ghost layers: one exchange feeds U (1) + curl (1) + curl
         *  (1) per sweep. MG coarsening keeps boxes >= 4 cells so periodic
         *  exchanges never wrap more than one period. */
        static constexpr int NGROW = 3;

        bool m_is_defined = false;
        bool m_nested = false;
        int m_mg_floor = 4;
        RT m_beta_max = 0.0;
        int m_nsweeps = 48;
        bool m_verbose = false;

        Ops* m_ops = nullptr;

        int m_nlevels = 1;
        amrex::Vector<amrex::Geometry> m_geom;

        /** Per level: l2 = sum_d 2/dx_d^2 (the diagonal of -del2). */
        std::vector<RT> m_l2;
        std::vector<amrex::GpuArray<amrex::Real, 3>> m_inv2dx;
        std::vector<amrex::GpuArray<amrex::Real, 3>> m_invdx2;

        /** Per-node frozen coefficients (per level):
         *  comp 0-2: bhat, comp 3: c_w, comp 4: a_e, comp 5: a_H,
         *  comp 6-8: w_c = Sigma_cc |B_tot| / rho_c (deposited-mass-matrix
         *  ion response through the Hall term; zero unless
         *  implicit_evolve.use_mass_matrices_pc = 1). */
        std::vector<std::unique_ptr<MF>> m_coef;

        /** Per level (all 3-component): W iterate (NGROW ghosts), RHS (0),
         *  residual (1, for restriction), and the two stage scratches. */
        std::vector<std::unique_ptr<MF>> m_xW;
        std::vector<std::unique_ptr<MF>> m_bW;
        std::vector<std::unique_ptr<MF>> m_rW;
        std::vector<std::unique_ptr<MF>> m_tmp;
        std::vector<std::unique_ptr<MF>> m_tmp2;

        /** Fine level: ghosted copy of the PC input b (2 ghosts: one
         *  exchange feeds the curl-curl RHS assembly). */
        std::unique_ptr<MF> m_bE;

        /** Inner-GMRES storage (nested mode): RHS, Arnoldi basis, scratch. */
        std::unique_ptr<MF> m_bK;
        std::unique_ptr<MF> m_z;
        std::unique_ptr<MF> m_w;
        std::vector<std::unique_ptr<MF>> m_kry;

        /** Physical-BC closure state: per index dimension and side, whether
         *  the face is a PEC(_Insulator) wall the closure must model; the
         *  index-dim -> vector-component map; and the EB mask (the nodal
         *  eb_update_E flag on the collocated grid; nullptr without EB). */
        bool m_has_pec_wall = false;
        /** RZ: the r-lo face is the symmetry axis (rmin == 0). */
        bool m_has_axis = false;
        bool m_wall_pec[AMREX_SPACEDIM][2] = {};
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
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::grid_type == ablastr::utils::enums::GridType::Collocated,
        "pc_hybrid_pic currently requires warpx.grid_type = collocated");

    // Physical-BC closure: detect PEC(_Insulator) faces so the PC's operator
    // reproduces the residual's boundary rows (identity wall-tangential-E
    // rows with dead columns; parity-mirrored ghosts for the delta-B /
    // delta-J stages). Other non-periodic BCs keep the previous
    // zero-extension closure (warned once).
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
                    ablastr::warn_manager::WMRecordWarning("HybridPICPC",
                        "pc_hybrid_pic: non-periodic, non-PEC field boundary "
                        "detected; the preconditioner keeps a homogeneous "
                        "zero-extension closure on that face (approximate "
                        "boundary rows).");
                }
            }
        }
    }

    const amrex::ParmParse pp("pc_hybrid_pic");
    std::string mode = "smoother";
    pp.query("mode", mode);
    pp.query("nsweeps", m_nsweeps);
    pp.query("verbose", m_verbose);
    pp.query("mg_floor", m_mg_floor);
    pp.query("sigma", m_sigma);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_sigma > 0.0 && m_sigma < 2.0,
        "pc_hybrid_pic.sigma must be in (0, 2)");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_mg_floor >= 4,
        "pc_hybrid_pic.mg_floor must be >= 4");
    if (mode == "nested") {
        m_nested = true;
        // The inner GMRES stops at an RHS-dependent iteration count, so the
        // preconditioner is nonlinear: classical (fixed-M) right-PC GMRES is
        // silently invalid. Enforce a flexible outer solver; the KSP type
        // itself is checked after KSPSetFromOptions in WarpX_PETSc.cpp.
        {
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
        }
        amrex::Print() << "pc_hybrid_pic: nested mode (MG + inner GMRES) -- the "
                       << "preconditioner is nonlinear; the outer Krylov solver must be "
                       << "flexible (newton.linear_solver = petsc_ksp, "
                       << "PETSC_OPTIONS=\"-ksp_type fgmres\")\n";
    } else {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(mode == "smoother",
            "pc_hybrid_pic.mode must be 'smoother' or 'nested'");
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_nsweeps > 0,
        "pc_hybrid_pic.nsweeps must be positive");

    // Storage on the solver-vector layout (nodal, single AMR level); in
    // nested mode the MG hierarchy coarsens that grid while every box stays
    // coarsenable, keeping the DistributionMapping so transfers are local.
    const amrex::MultiFab* E0 = a_U.getArrayVec()[0][0];
    const amrex::DistributionMapping& dm = E0->DistributionMap();

    amrex::Vector<amrex::BoxArray> ba_node;
    ba_node.push_back(E0->boxArray());
    m_geom.clear();
    m_geom.push_back(warpx.Geom(0));

    if (m_nested) {
        amrex::BoxArray ba_cell = amrex::convert(ba_node[0], amrex::IntVect::TheCellVector());
        while (true) {
            const amrex::Geometry& fg = m_geom.back();
            // min box width 4: periodic exchanges of NGROW ghosts must not
            // wrap more than one period
            if (!fg.Domain().coarsenable(2, m_mg_floor)
                || !ba_cell.coarsenable(2, m_mg_floor)) { break; }
            ba_cell.coarsen(2);
            ba_node.push_back(amrex::convert(ba_cell, amrex::IntVect::TheNodeVector()));
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
    m_inv2dx.assign(m_nlevels, {0.0, 0.0, 0.0});
    m_invdx2.assign(m_nlevels, {0.0, 0.0, 0.0});
    for (int lev = 0; lev < m_nlevels; lev++) {
        const auto* dx = m_geom[lev].CellSize();
        for (int d = 0; d < ndirs; d++) {
            const amrex::Real h = dx[d];
            m_inv2dx[lev][dirs[d]] = 0.5_rt / h;
            m_invdx2[lev][dirs[d]] = 1.0_rt / (h*h);
#if defined(WARPX_DIM_RZ)
            // the RZ curl-curl is the compact down/up composition, whose
            // Nyquist symbol reaches 4/h^2 per direction (the Cartesian
            // wide central composition peaks at half that): size the
            // defect diagonal accordingly so the sigma < 2 damping bound
            // keeps its Cartesian semantics
            m_l2[lev] += 4.0_rt / (h*h);
#else
            m_l2[lev] += 2.0_rt / (h*h);
#endif
        }
    }

    m_coef.resize(m_nlevels);
    m_xW.resize(m_nlevels);
    m_bW.resize(m_nlevels);
    m_rW.resize(m_nlevels);
    m_tmp.resize(m_nlevels);
    m_tmp2.resize(m_nlevels);

    for (int lev = 0; lev < m_nlevels; lev++) {
        const amrex::BoxArray& ba = ba_node[lev];
        // 2 coefficient ghost layers: the fused U-stage evaluates U in a
        // 2-ring around the valid region; internal ghosts are exchanged in
        // Update(), physical ghosts stay zero (U vanishes there, matching
        // the homogeneous-Dirichlet closure of the unfused operator)
        m_coef[lev] = std::make_unique<MF>(ba, dm, 9, 2);
        m_coef[lev]->setVal(0.0);
        m_xW[lev]   = std::make_unique<MF>(ba, dm, 3, NGROW);
        m_rW[lev]   = std::make_unique<MF>(ba, dm, 3, 1);
        m_tmp[lev]  = std::make_unique<MF>(ba, dm, 3, NGROW);
        m_tmp2[lev] = std::make_unique<MF>(ba, dm, 3, NGROW);
        m_bW[lev]   = std::make_unique<MF>(ba, dm, 3, 0);
        // Domain (non-periodic) ghosts stay zero: homogeneous Dirichlet
        // closure on every level (exchanges fill periodic/internal only)
        m_xW[lev]->setVal(0.0);
        m_rW[lev]->setVal(0.0);
        m_tmp[lev]->setVal(0.0);
        m_tmp2[lev]->setVal(0.0);
        m_bW[lev]->setVal(0.0);
    }
    m_bE = std::make_unique<MF>(ba_node[0], dm, 3, 2);
    m_bE->setVal(0.0);

    if (m_nested) {
        m_kry.resize(INNER_MAX + 1);
        m_bK = std::make_unique<MF>(ba_node[0], dm, 3, 0);
        m_w  = std::make_unique<MF>(ba_node[0], dm, 3, 0);
        m_z  = std::make_unique<MF>(ba_node[0], dm, 3, NGROW);
        m_bK->setVal(0.0);
        m_w->setVal(0.0);
        m_z->setVal(0.0);
        for (int j = 0; j <= INNER_MAX; j++) {
            m_kry[j] = std::make_unique<MF>(ba_node[0], dm, 3, 0);
            m_kry[j]->setVal(0.0);
        }
    }

    m_is_defined = true;
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::printParameters () const
{
    using namespace amrex;
    Print() << "Hybrid-PIC EMHD preconditioner (pc_hybrid_pic, reduced W-system):\n";
    if (m_nested) {
        Print() << "  mode: nested -- MG V(" << NU_PRE << "," << NU_POST << ") ("
                << NU_BOTTOM << " bottom sweeps, coarsest "
                << m_geom.back().Domain().length(0) << ") inside GMRES(<= "
                << INNER_MAX << " its, W-space rtol " << INNER_RTOL
                << "/max(1, beta_max) targeting E-space " << INNER_RTOL
                << "); requires a flexible outer solver\n";
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
    Print() << "  damping: " << m_sigma << "\n";
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
    const amrex::Periodicity& period0 = m_geom[0].periodicity();
    for (int n = 0; n < 3; n++) {
        amrex::MultiFab::Copy(*m_bE, *B_fp[n], 0, n, 1, 0);
        if (add_external) {
            amrex::MultiFab::Add(*m_bE, *B_ext[n], 0, n, 1, 0);
        }
    }
    m_bE->FillBoundary(period0);
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
        // on the collocated grid all J components share the coef's nodal
        // layout; anything else means an unsupported staggering slipped in
        AMREX_ALWAYS_ASSERT(Sxx_mf->boxArray() == m_coef[0]->boxArray()
            && Sxx_mf->DistributionMap() == m_coef[0]->DistributionMap());
    }

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
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real bx_v = Bx(i,j,k);
                Real by_v = By(i,j,k);
                Real bz_v = Bz(i,j,k);
                if (add_external) {
                    bx_v += Bxe(i,j,k);
                    by_v += Bye(i,j,k);
                    bz_v += Bze(i,j,k);
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
            });
    }

    // Peak whistler CFL: calibrates the inner-GMRES W-space tolerance so
    // it targets E-space accuracy (r_E = -U r_W; see the class docs)
    if (m_nested) {
        m_beta_max = static_cast<RT>(m_l2[0] * m_coef[0]->max(3, 0));
    }

    // Coarse operators by rediscretization: inject the coefficients down the
    // hierarchy (l2 is recomputed per level in Define()) and exchange their
    // internal ghosts for the fused U-stage.
    m_coef[0]->FillBoundary(period0);
    ApplyDomainClosure(*m_coef[0], 0, ClosureClass::Coef);
    for (int lev = 1; lev < m_nlevels; lev++) {
        RestrictCoef(lev);
        m_coef[lev]->FillBoundary(m_geom[lev].periodicity());
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
        amrex::ParallelFor(bx, 9,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
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
                                       CurlKind a_kind ) const
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
        Array4<const Real> const& av = a.const_array(mfi);
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
    amrex::ignore_unused(a_kind);
    const auto inv2dx = m_inv2dx[a_lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(a_grow);
        Array4<const Real> const& av = a.const_array(mfi);
        Array4<Real> const& cv = out.array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
#if defined(WARPX_DIM_3D)
                const Real dyaz = (av(i,j+1,k,2) - av(i,j-1,k,2)) * inv2dx[1];
                const Real dzay = (av(i,j,k+1,1) - av(i,j,k-1,1)) * inv2dx[2];
                const Real dzax = (av(i,j,k+1,0) - av(i,j,k-1,0)) * inv2dx[2];
                const Real dxaz = (av(i+1,j,k,2) - av(i-1,j,k,2)) * inv2dx[0];
                const Real dxay = (av(i+1,j,k,1) - av(i-1,j,k,1)) * inv2dx[0];
                const Real dyax = (av(i,j+1,k,0) - av(i,j-1,k,0)) * inv2dx[1];
                cv(i,j,k,0) = dyaz - dzay;
                cv(i,j,k,1) = dzax - dxaz;
                cv(i,j,k,2) = dxay - dyax;
#elif defined(WARPX_DIM_XZ)
                // (i,j) = (x,z); d/dy = 0
                const Real dzay = (av(i,j+1,k,1) - av(i,j-1,k,1)) * inv2dx[2];
                const Real dzax = (av(i,j+1,k,0) - av(i,j-1,k,0)) * inv2dx[2];
                const Real dxaz = (av(i+1,j,k,2) - av(i-1,j,k,2)) * inv2dx[0];
                const Real dxay = (av(i+1,j,k,1) - av(i-1,j,k,1)) * inv2dx[0];
                cv(i,j,k,0) = -dzay;
                cv(i,j,k,1) = dzax - dxaz;
                cv(i,j,k,2) = dxay;
#else // WARPX_DIM_1D_Z
                // i = z; d/dx = d/dy = 0
                const Real dzay = (av(i+1,j,k,1) - av(i-1,j,k,1)) * inv2dx[2];
                const Real dzax = (av(i+1,j,k,0) - av(i-1,j,k,0)) * inv2dx[2];
                cv(i,j,k,0) = -dzay;
                cv(i,j,k,1) = dzax;
                cv(i,j,k,2) = 0.0_rt;
#endif
            });
    }
#endif // WARPX_DIM_RZ
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyUValid ( const MF& a, MF& out, int a_lev ) const
{
    using namespace amrex;

    const auto invdx2 = m_invdx2[a_lev];
#if defined(WARPX_DIM_RZ)
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
        const Box& bx = mfi.tilebox();
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
                ov(i,j,k,0) = cw*(ay*bhz - az*bhy) + ae*ax - ah*lx;
                ov(i,j,k,1) = cw*(az*bhx - ax*bhz) + ae*ay - ah*ly;
                ov(i,j,k,2) = cw*(ax*bhy - ay*bhx) + ae*az - ah*lz;
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::StageKU ( const MF& a, int a_lev )
{
    using namespace amrex;

    const auto invdx2 = m_invdx2[a_lev];
    const bool ion = m_have_mm;
#if defined(WARPX_DIM_RZ)
    const auto* dxg = m_geom[a_lev].CellSize();
    const Real rz_dr = dxg[0];
    const Real rz_inv_dr = 1.0_rt/dxg[0];
    const Real rz_inv_dz = 1.0_rt/dxg[1];
    const Real rz_rmin = static_cast<Real>(m_geom[a_lev].ProbLo(0));
#endif

    // tmp = U a on growntilebox(2): U reads a's first ghost ring beyond
    // that, so a must carry NGROW = 3 filled ghosts
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_tmp[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox(2);
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
                Real ux = cw*(ay*bhz - az*bhy) + ae*ax - ah*lx;
                Real uy = cw*(az*bhx - ax*bhz) + ae*ay - ah*ly;
                Real uz = cw*(ax*bhy - ay*bhx) + ae*az - ah*lz;
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

    // tmp2 = curl tmp on growntilebox(1)
    ComputeCurl(*m_tmp[a_lev], *m_tmp2[a_lev], a_lev, 1, CurlKind::EtoB);
    ApplyDomainClosure(*m_tmp2[a_lev], a_lev, ClosureClass::BField);
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
    const auto inv2dx = m_inv2dx[a_lev];
    const amrex::Periodicity& period = m_geom[a_lev].periodicity();

    for (int sweep = 0; sweep < a_nsweeps; sweep++) {

        const bool r_is_b = (sweep == 0 && a_x_is_zero);
        if (!r_is_b) {
            // one exchange feeds the whole U -> curl -> curl composition
            m_xW[a_lev]->FillBoundary(period);
            // W is delta-J-like: parity ghosts, no wall write
            ApplyDomainClosure(*m_xW[a_lev], a_lev, ClosureClass::JField);
            StageKU(*m_xW[a_lev], a_lev);
        }

        {
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
                            const Real ccx = (t2(i,j+1,k,2) - t2(i,j-1,k,2))*inv2dx[1]
                                           - (t2(i,j,k+1,1) - t2(i,j,k-1,1))*inv2dx[2];
                            const Real ccy = (t2(i,j,k+1,0) - t2(i,j,k-1,0))*inv2dx[2]
                                           - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                            const Real ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0]
                                           - (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[1];
#elif defined(WARPX_DIM_XZ)
                            const Real ccx = -(t2(i,j+1,k,1) - t2(i,j-1,k,1))*inv2dx[2];
                            const Real ccy = (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[2]
                                           - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                            const Real ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0];
#else
                            const Real ccx = -(t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[2];
                            const Real ccy = (t2(i+1,j,k,0) - t2(i-1,j,k,0))*inv2dx[2];
                            const Real ccz = 0.0_rt;
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
                        // isotropic (rotation) part, |S| ~ Omega_ci theta*dt
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

    const amrex::Periodicity& period = m_geom[a_lev].periodicity();
    const auto inv2dx = m_inv2dx[a_lev];
    const auto invdx2 = m_invdx2[a_lev];

    m_xW[a_lev]->FillBoundary(period);
    ApplyDomainClosure(*m_xW[a_lev], a_lev, ClosureClass::JField);
    StageKU(*m_xW[a_lev], a_lev);

    {
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
                    const Real ccx = (t2(i,j+1,k,2) - t2(i,j-1,k,2))*inv2dx[1]
                                   - (t2(i,j,k+1,1) - t2(i,j,k-1,1))*inv2dx[2];
                    const Real ccy = (t2(i,j,k+1,0) - t2(i,j,k-1,0))*inv2dx[2]
                                   - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                    const Real ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0]
                                   - (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[1];
#elif defined(WARPX_DIM_XZ)
                    const Real ccx = -(t2(i,j+1,k,1) - t2(i,j-1,k,1))*inv2dx[2];
                    const Real ccy = (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[2]
                                   - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                    const Real ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0];
#else
                    const Real ccx = -(t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[2];
                    const Real ccy = (t2(i+1,j,k,0) - t2(i-1,j,k,0))*inv2dx[2];
                    const Real ccz = 0.0_rt;
#endif
                    rW(i,j,k,0) = bW(i,j,k,0) - xW(i,j,k,0) - ccx;
                    rW(i,j,k,1) = bW(i,j,k,1) - xW(i,j,k,1) - ccy;
                    rW(i,j,k,2) = bW(i,j,k,2) - xW(i,j,k,2) - ccz;
                });
        }
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyReducedOp ( const MF& a, MF& out )
{
    using namespace amrex;
#if defined(WARPX_DIM_RZ)
    const auto* dxg = m_geom[0].CellSize();
    const Real rz_dr = dxg[0];
    const Real rz_inv_dr = 1.0_rt/dxg[0];
    const Real rz_inv_dz = 1.0_rt/dxg[1];
    const Real rz_rmin = static_cast<Real>(m_geom[0].ProbLo(0));
#endif

    // a must be ghost-filled (NGROW layers) by the caller
    StageKU(a, 0);
    const auto inv2dx = m_inv2dx[0];
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
                const Real ccx = (t2(i,j+1,k,2) - t2(i,j-1,k,2))*inv2dx[1]
                               - (t2(i,j,k+1,1) - t2(i,j,k-1,1))*inv2dx[2];
                const Real ccy = (t2(i,j,k+1,0) - t2(i,j,k-1,0))*inv2dx[2]
                               - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                const Real ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0]
                               - (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[1];
#elif defined(WARPX_DIM_XZ)
                const Real ccx = -(t2(i,j+1,k,1) - t2(i,j-1,k,1))*inv2dx[2];
                const Real ccy = (t2(i,j+1,k,0) - t2(i,j-1,k,0))*inv2dx[2]
                               - (t2(i+1,j,k,2) - t2(i-1,j,k,2))*inv2dx[0];
                const Real ccz = (t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[0];
#else
                const Real ccx = -(t2(i+1,j,k,1) - t2(i-1,j,k,1))*inv2dx[2];
                const Real ccy = (t2(i+1,j,k,0) - t2(i-1,j,k,0))*inv2dx[2];
                const Real ccz = 0.0_rt;
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_bW[a_lev+1], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<const Real> const& fr = m_rW[a_lev]->const_array(mfi);
        Array4<Real> const& cb = m_bW[a_lev+1]->array(mfi);
        amrex::ParallelFor(bx, 3,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                // nodal full weighting; fine ghosts are filled (periodic/
                // internal) or zero (physical: homogeneous Dirichlet)
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
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ProlongLevel ( int a_lev )
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::ProlongLevel()");

    const amrex::Periodicity& cperiod = m_geom[a_lev+1].periodicity();
    m_xW[a_lev+1]->FillBoundary(cperiod);
    ApplyDomainClosure(*m_xW[a_lev+1], a_lev+1, ClosureClass::JField);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_xW[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<const Real> const& xc = m_xW[a_lev+1]->const_array(mfi);
        Array4<Real> const& xf = m_xW[a_lev]->array(mfi);
        amrex::ParallelFor(bx, 3,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                // (tri)linear interpolation on nodes (level-0 index space
                // starts at 0, so i/2 floors correctly)
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
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::Vcycle ( int a_lev, bool a_x_is_zero )
{
    if (a_lev == m_nlevels - 1) {
        SmoothLevel(a_lev, (m_nlevels == 1) ? m_nsweeps : NU_BOTTOM, a_x_is_zero);
        return;
    }

    SmoothLevel(a_lev, NU_PRE, a_x_is_zero);

    ResidualLevel(a_lev);
    m_rW[a_lev]->FillBoundary(m_geom[a_lev].periodicity());
    ApplyDomainClosure(*m_rW[a_lev], a_lev, ClosureClass::JField);
    RestrictLevel(a_lev);
    m_xW[a_lev+1]->setVal(0.0);
    Vcycle(a_lev + 1, true);
    ProlongLevel(a_lev);

    SmoothLevel(a_lev, NU_POST, false);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyMToW ( const MF& a_in )
{
    using namespace amrex;

    amrex::MultiFab::Copy(*m_bW[0], a_in, 0, 0, 3, 0);
    m_xW[0]->setVal(0.0);
    Vcycle(0, true);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::RotateValid ( MF& a ) const
{
    using namespace amrex;
    // a <- (I+S)^-1 a pointwise (deposited-mass-matrix ion response,
    // see class docs); valid region only
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
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
    RotateValid(*m_bE);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ApplyDomainClosure ( MF& a, int a_lev,
                                              ClosureClass a_cls ) const
{
    using namespace amrex;

    if (m_has_pec_wall) {
        const Box dom = amrex::convert(m_geom[a_lev].Domain(),
                                       a.boxArray().ixType());
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            for (int side = 0; side < 2; side++) {
                if (!m_wall_pec[d][side]) { continue; }
                const int w = (side == 0) ? dom.smallEnd(d) : dom.bigEnd(d);
                const int nc = m_dim2comp[d];
                const bool lo_side = (side == 0);
                const int cls = static_cast<int>(a_cls);
                const int dd = d;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.growntilebox();
                    Array4<Real> const& v = a.array(mfi);
                    amrex::ParallelFor(bx, a.nComp(),
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            const int id = (dd == 0) ? i : ((dd == 1) ? j : k);
                            const bool outside =
                                lo_side ? (id < w) : (id > w);
                            // Coef: comps 0-2 carry bhat (vector), 3-8 scalar
                            const bool scalar =
                                (cls == 4/*Coef*/) && (n >= 3);
                            const bool is_norm = !scalar && (n == nc);
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
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox();
            Array4<Real> const& v = a.array(mfi);
            amrex::ParallelFor(bx, a.nComp(),
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    if (i > w) { return; }
                    const bool scalar = (cls == 4/*Coef*/) && (n >= 3);
                    if (i == w) {
                        if (cls == 0/*EField*/ && n == 1) {
                            v(i,j,k,n) = 0.0_rt;
                        }
                        return;
                    }
                    const int im = 2*w - i;
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
    if (a_lev == 0 && m_eb_flag != nullptr) {
        amrex::IntVect ngm = a.nGrowVect();
        ngm.min(m_eb_flag->nGrowVect());
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(ngm);
            Array4<Real> const& v = a.array(mfi);
            Array4<const int> const& fl = m_eb_flag->const_array(mfi);
            amrex::ParallelFor(bx, a.nComp(),
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    if (fl(i,j,k) == 0) { v(i,j,k,n) = 0.0_rt; }
                });
        }
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::OverrideIdentityRows ( T& a_out, const T& a_in ) const
{
    using namespace amrex;

    if (m_has_pec_wall) {
        const Box dom = amrex::convert(m_geom[0].Domain(),
            a_out.getArrayVec()[0][0]->boxArray().ixType());
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            for (int side = 0; side < 2; side++) {
                if (!m_wall_pec[d][side]) { continue; }
                const int w = (side == 0) ? dom.smallEnd(d) : dom.bigEnd(d);
                Box plane = dom;
                plane.setSmall(d, w);
                plane.setBig(d, w);
                const int nc = m_dim2comp[d];
                for (int n = 0; n < 3; n++) {
                    if (n == nc) { continue; }  // tangential rows only
                    MF& o = *a_out.getArrayVec()[0][n];
                    const MF& s = *a_in.getArrayVec()[0][n];
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
typename HybridPICPC<T,Ops>::RT HybridPICPC<T,Ops>::DotW ( const MF& a, const MF& b ) const
{
    using warpx::fields::FieldType;
    using ablastr::fields::Direction;

    auto& warpx = WarpX::GetInstance();
    // On the collocated grid all three components share the nodal Efield_fp
    // staggering: one dot mask serves all comps (no double counting of
    // box-shared nodes)
    const amrex::iMultiFab* mask =
        warpx.getFieldDotMaskPointer(FieldType::Efield_fp, 0, Direction{0});
    amrex::Real result = amrex::MultiFab::Dot(*mask, a, 0, b, 0, 3, 0, true);
    amrex::ParallelAllReduce::Sum(result, amrex::ParallelContext::CommunicatorSub());
    return static_cast<RT>(result);
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::InnerGMRES ()
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::InnerGMRES()");

    const int m = INNER_MAX;
    // W-space tolerance scaled so the E-space residual (r_E = -U r_W, with
    // ||U|| up to beta_max at the grid scale) meets INNER_RTOL
    const RT rtol = static_cast<RT>(INNER_RTOL)
                  / amrex::max(RT(1.0), m_beta_max);
    const amrex::Periodicity& period = m_geom[0].periodicity();

    amrex::MultiFab::Copy(*m_kry[0], *m_bK, 0, 0, 3, 0);
    m_xW[0]->setVal(0.0);
    const RT beta = std::sqrt(DotW(*m_kry[0], *m_kry[0]));
    if (beta <= RT(0.0)) { return; }
    m_kry[0]->mult(RT(1.0)/beta, 0, 3, 0);

    std::vector<std::vector<RT>> H(m + 1, std::vector<RT>(m, RT(0.0)));
    std::vector<RT> cs(m, RT(0.0));
    std::vector<RT> sn(m, RT(0.0));
    std::vector<RT> g(m + 1, RT(0.0));
    g[0] = beta;

    int its = 0;
    for (int j = 0; j < m; j++) {
        its = j + 1;

        // z = M v_j (one V-cycle; result in m_xW[0]), w = A z
        ApplyMToW(*m_kry[j]);
        amrex::MultiFab::Copy(*m_z, *m_xW[0], 0, 0, 3, 0);
        m_z->FillBoundary(period);
        ApplyDomainClosure(*m_z, 0, ClosureClass::JField);
        ApplyReducedOp(*m_z, *m_w);

        for (int i = 0; i <= j; i++) {
            H[i][j] = DotW(*m_w, *m_kry[i]);
            amrex::MultiFab::Saxpy(*m_w, -H[i][j], *m_kry[i], 0, 0, 3, 0);
        }
        const RT hj1 = std::sqrt(DotW(*m_w, *m_w));

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

        amrex::MultiFab::Copy(*m_kry[j+1], *m_w, 0, 0, 3, 0);
        m_kry[j+1]->mult(RT(1.0)/hj1, 0, 3, 0);
    }

    std::vector<RT> y(its, RT(0.0));
    for (int i = its - 1; i >= 0; i--) {
        RT s = g[i];
        for (int k = i + 1; k < its; k++) { s -= H[i][k]*y[k]; }
        y[i] = s / H[i][i];
    }

    // x = M (V y): the V-cycle is linear, so apply it once to the summed
    // basis (classical right-preconditioned GMRES update)
    m_z->setVal(0.0);
    for (int i = 0; i < its; i++) {
        amrex::MultiFab::Saxpy(*m_z, y[i], *m_kry[i], 0, 0, 3, 0);
    }
    ApplyMToW(*m_z);   // final W left in m_xW[0]

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

    const amrex::Periodicity& period = m_geom[0].periodicity();

    // bW = curl curl b (exact operators, one exchange): the curl-free
    // content of b never enters the W-solve and passes through the
    // reconstruction exactly.
    for (int n = 0; n < 3; n++) {
        amrex::MultiFab::Copy(*m_bE, *a_b.getArrayVec()[0][n], 0, n, 1, 0);
    }
    if (m_have_mm) {
        // the W-system RHS is K (I+S)^-1 b: rotate b in place first (the
        // reconstruction below re-reads m_bE, which must stay the ROTATED
        // b so that x = (I+S)^-1 b - (I+S)^-1 U W)
        RotateBInPlace();
    }
    m_bE->FillBoundary(period);
    // PEC dead tangential-E columns (wall values zeroed; the identity rows
    // are restored on the output below) + parity ghosts; EB nodes zeroed
    ApplyDomainClosure(*m_bE, 0, ClosureClass::EField);
    m_xW[0]->setVal(0.0);
    ComputeCurl(*m_bE, *m_tmp[0], 0, 1, CurlKind::EtoB);   // curl b on grow(1)
    ApplyDomainClosure(*m_tmp[0], 0, ClosureClass::BField);
    MF& bdst = m_nested ? *m_bK : *m_bW[0];
    ComputeCurl(*m_tmp[0], bdst, 0, 0, CurlKind::BtoJ);    // curl curl b (valid)

    if (m_nested) {
        InnerGMRES();
    } else {
        SmoothLevel(0, m_nsweeps, true);
    }

    // Reconstruct x = (I+S)^-1 (b - U xW) = bE_rot - (I+S)^-1 U xW
    // (m_bE already holds the rotated b when the ion-response leg is on)
    m_xW[0]->FillBoundary(period);
    ApplyDomainClosure(*m_xW[0], 0, ClosureClass::JField);
    ApplyUValid(*m_xW[0], *m_tmp[0], 0);
    if (m_have_mm) { RotateValid(*m_tmp[0]); }
    // the residual BCs the Ohm output: wall-tangential/EB contributions off
    ApplyDomainClosure(*m_tmp[0], 0, ClosureClass::EField);
    for (int n = 0; n < 3; n++) {
        amrex::MultiFab::LinComb(*a_x.getArrayVec()[0][n],
                                 1.0_rt, *m_bE, n,
                                 -1.0_rt, *m_tmp[0], n,
                                 0, 1, 0);
    }
    // exact identity rows of the true Jacobian: x = b on wall-tangential E
    // at PEC faces and on EB-covered nodes
    OverrideIdentityRows(a_x, a_b);
}

#endif
