/* Copyright 2026 Prabhat Kumar
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_HYBRID_PIC_PC_H_
#define WARPX_HYBRID_PIC_PC_H_

#include "Preconditioner.H"

#include "Fields.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/utils/Enums.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>

/**
 * \brief Physics-based preconditioner for the theta-implicit hybrid-PIC solver
 *  (electron-MHD block), following L. Chacon, JCP 530 (2025) 113910.
 *
 *  The stiff part of the hybrid Ohm's-law Jacobian is the electron (EMHD)
 *  block: with particles, density and electron pressure frozen, the only
 *  dependence of E_Ohm on the Newton iterate E is through Faraday's law,
 *  dB = -theta*dt curl(dE), giving
 *
 *    J_pc dE = dE + theta*dt [ (curl curl dE) x B0 / (mu0 e n)     (whistler)
 *                              + (eta/mu0) curl curl dE            (resistive)
 *                              - (eta_h/mu0) del2 (curl curl dE) ] (hyper-res.)
 *
 *  Following the paper, the 4th-order system is split into two 2nd-order ones
 *  by introducing the auxiliary unknown W = curl curl E (proportional to the
 *  current-density response):
 *
 *    Row 1:  dE + theta*dt [ cw (dW x bhat) + (eta/mu0) dW - (eta_h/mu0) del2 dW ] = b
 *    Row 2:  dW - curl curl dE = 0,       cw = |B0| / (mu0 e n)  [times theta*dt here]
 *
 *  and solved approximately with damped block-Jacobi sweeps. The per-node
 *  6x6 block diagonal is built with the defect-correction approximation
 *  curl curl ~ -del2 (diagonal l2 = sum_d 2/dx_d^2), which has the analytic
 *  inverse via P = alpha I + beta (I x bhat):
 *
 *    P^-1 v = [ alpha v - beta (v x bhat) + (beta^2/alpha)(v . bhat) bhat ] / (alpha^2 + beta^2)
 *    alpha  = 1 + l2 * aU,   beta = l2 * cw,   aU = theta*dt*(eta + eta_h*l2)/mu0
 *
 *  All V-cycle-style residuals use the exact central-difference operators, so
 *  the approximation only affects the iteration, not its fixed point.
 *
 *  This version supports the collocated (nodal) grid in Cartesian 1D/2D/3D,
 *  where the per-node block is exact (E, B, rho co-located). Stability of the
 *  damped Jacobi iteration requires sigma < 2/(max l2 dx^2 / 2 ... ): at the
 *  Nyquist mode lambda = 1 - sigma*kappa2 with kappa2 <= AMREX_SPACEDIM, so
 *  sigma <= 0.6 is used by default in 3D (the paper's 0.7 is a 2D value).
 *
 *  Input parameters (prefix pc_hybrid_pic.):
 *    nsweeps  (int,  default 6)    number of damped block-Jacobi sweeps
 *    damping  (Real, default 0.55) Jacobi damping parameter sigma
 *    verbose  (bool, default false)
 */
template <class T, class Ops>
class HybridPICPC : public Preconditioner<T,Ops>
{
    public:

        using RT = typename T::value_type;

        HybridPICPC () = default;
        ~HybridPICPC () override = default;

        HybridPICPC(const HybridPICPC&) = delete;
        HybridPICPC& operator=(const HybridPICPC&) = delete;
        HybridPICPC(HybridPICPC&&) noexcept = delete;
        HybridPICPC& operator=(HybridPICPC&&) noexcept = delete;

        void Define (const T&, Ops*) override;

        /** Recompute the frozen coefficients (bhat, cw, eta terms) from the
         *  current linearization state (B^{n+theta} in Bfield_fp, deposited
         *  rho in rho_fp). Called once per Newton iteration. */
        void Update (const T& a_U) override;

        /** Approximately solve J_pc x = b with damped block-Jacobi sweeps on
         *  the split (dE, dW) system. */
        void Apply (T& a_x, const T& a_b) override;

        void printParameters () const override;

        [[nodiscard]] inline bool IsDefined () const override { return m_is_defined; }

        /** Central-difference curl of the ghosted fields a into out (valid region).
         * Public because the CUDA compiler requires the enclosing function of an
         * extended device lambda to have public access. */
        void ComputeCurl (const std::array<std::unique_ptr<amrex::MultiFab>, 3>& a,
                          std::array<std::unique_ptr<amrex::MultiFab>, 3>& out) const;

    private:

        bool m_is_defined = false;

        int m_nsweeps = 6;
        RT m_sigma = 0.55;
        bool m_verbose = false;

        Ops* m_ops = nullptr;

        /** l2 = sum_d 2/dx_d^2, the (positive) diagonal of -del2 */
        RT m_l2 = 0.0;
        amrex::GpuArray<amrex::Real, 3> m_inv2dx = {0.0, 0.0, 0.0};
        amrex::GpuArray<amrex::Real, 3> m_invdx2 = {0.0, 0.0, 0.0};

        /** Per-node frozen coefficients:
         *  comp 0-2: bhat (unit B, zero vector where |B| = 0)
         *  comp 3:   cw    = theta*dt*|B|/(mu0*rho_e)   (whistler)
         *  comp 4:   a_eta = theta*dt*eta/mu0           (resistive)
         *  comp 5:   a_H   = theta*dt*eta_h/mu0         (hyper-resistive) */
        std::unique_ptr<amrex::MultiFab> m_coef;

        /** Scratch: iterate (ghosted for stencils) and residuals */
        std::array<std::unique_ptr<amrex::MultiFab>, 3> m_dE;
        std::array<std::unique_ptr<amrex::MultiFab>, 3> m_dW;
        std::array<std::unique_ptr<amrex::MultiFab>, 3> m_tmp;
        std::array<std::unique_ptr<amrex::MultiFab>, 3> m_r1;
        std::array<std::unique_ptr<amrex::MultiFab>, 3> m_r2;

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
#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    WARPX_ABORT_WITH_MESSAGE("pc_hybrid_pic is not implemented for curvilinear geometries");
#endif
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::grid_type == ablastr::utils::enums::GridType::Collocated,
        "pc_hybrid_pic currently requires warpx.grid_type = collocated");

    const amrex::ParmParse pp("pc_hybrid_pic");
    pp.query("nsweeps", m_nsweeps);
    pp.query("damping", m_sigma);
    pp.query("verbose", m_verbose);

    // 3D Nyquist stability bound of the damped Jacobi iteration:
    // lambda = 1 - sigma*kappa2 with kappa2 up to AMREX_SPACEDIM.
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_sigma > 0.0 && m_sigma * AMREX_SPACEDIM < 2.0,
        "pc_hybrid_pic.damping must satisfy 0 < sigma < 2/ndims for a stable Jacobi iteration");

    const auto* dx = warpx.Geom(0).CellSize();
    m_l2 = 0.0;
#if defined(WARPX_DIM_3D)
    const int dirs[3] = {0, 1, 2};
    constexpr int ndirs = 3;
#elif defined(WARPX_DIM_XZ)
    const int dirs[2] = {0, 2};
    constexpr int ndirs = 2;
#else // WARPX_DIM_1D_Z
    const int dirs[1] = {2};
    constexpr int ndirs = 1;
#endif
    for (int d = 0; d < ndirs; d++) {
        // AMReX dimension index of WarpX direction dirs[d]
        const amrex::Real h = dx[d];
        m_inv2dx[dirs[d]] = 0.5_rt / h;
        m_invdx2[dirs[d]] = 1.0_rt / (h*h);
        m_l2 += 2.0_rt / (h*h);
    }

    // All scratch on the solver-vector layout (nodal, single level)
    const amrex::MultiFab* E0 = a_U.getArrayVec()[0][0];
    const amrex::BoxArray& ba = E0->boxArray();
    const amrex::DistributionMapping& dm = E0->DistributionMap();
    const amrex::IntVect ng(1);

    m_coef = std::make_unique<MultiFab>(ba, dm, 6, 0);
    for (int n = 0; n < 3; n++) {
        m_dE[n]  = std::make_unique<MultiFab>(ba, dm, 1, ng);
        m_dW[n]  = std::make_unique<MultiFab>(ba, dm, 1, ng);
        m_tmp[n] = std::make_unique<MultiFab>(ba, dm, 1, ng);
        m_r1[n]  = std::make_unique<MultiFab>(ba, dm, 1, 0);
        m_r2[n]  = std::make_unique<MultiFab>(ba, dm, 1, 0);
        // Domain (non-periodic) ghosts stay zero: homogeneous Dirichlet in the PC
        m_dE[n]->setVal(0.0);
        m_dW[n]->setVal(0.0);
        m_tmp[n]->setVal(0.0);
    }

    m_is_defined = true;
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::printParameters () const
{
    using namespace amrex;
    Print() << "Hybrid-PIC EMHD preconditioner (pc_hybrid_pic):\n";
    Print() << "  nsweeps: " << m_nsweeps << "\n";
    Print() << "  damping: " << m_sigma << "\n";
    Print() << "  diag l2: " << m_l2 << "\n";
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
    const RT one_over_mu0 = 1.0_rt / PhysConst::mu0;
    const auto eta = hybrid->m_eta;
    const auto eta_h = hybrid->m_eta_h;
    const bool include_hyper = hybrid->m_include_hyper_resistivity_term;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_coef, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& coef = m_coef->array(mfi);
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
                const Real rho_e = amrex::max(r(i,j,k,0), rho_floor);
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
                // Linearize the resistivity at |J| = 0 (constant-eta production
                // runs are exact; J-dependent eta is only approximated)
                coef(i,j,k,4) = theta_dt * eta(rho_e, 0.0_rt, t_now) * one_over_mu0;
                coef(i,j,k,5) = include_hyper ?
                    theta_dt * eta_h(rho_e, Bmag) * one_over_mu0 : 0.0_rt;
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::ComputeCurl (
    const std::array<std::unique_ptr<amrex::MultiFab>, 3>& a,
    std::array<std::unique_ptr<amrex::MultiFab>, 3>& out ) const
{
    using namespace amrex;

    const auto inv2dx = m_inv2dx;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*out[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<const Real> const& ax = a[0]->const_array(mfi);
        Array4<const Real> const& ay = a[1]->const_array(mfi);
        Array4<const Real> const& az = a[2]->const_array(mfi);
        Array4<Real> const& cx = out[0]->array(mfi);
        Array4<Real> const& cy = out[1]->array(mfi);
        Array4<Real> const& cz = out[2]->array(mfi);
        amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
#if defined(WARPX_DIM_3D)
                const Real dyaz = (az(i,j+1,k) - az(i,j-1,k)) * inv2dx[1];
                const Real dzay = (ay(i,j,k+1) - ay(i,j,k-1)) * inv2dx[2];
                const Real dzax = (ax(i,j,k+1) - ax(i,j,k-1)) * inv2dx[2];
                const Real dxaz = (az(i+1,j,k) - az(i-1,j,k)) * inv2dx[0];
                const Real dxay = (ay(i+1,j,k) - ay(i-1,j,k)) * inv2dx[0];
                const Real dyax = (ax(i,j+1,k) - ax(i,j-1,k)) * inv2dx[1];
                cx(i,j,k) = dyaz - dzay;
                cy(i,j,k) = dzax - dxaz;
                cz(i,j,k) = dxay - dyax;
#elif defined(WARPX_DIM_XZ)
                // (i,j) = (x,z); d/dy = 0
                const Real dzay = (ay(i,j+1,k) - ay(i,j-1,k)) * inv2dx[2];
                const Real dzax = (ax(i,j+1,k) - ax(i,j-1,k)) * inv2dx[2];
                const Real dxaz = (az(i+1,j,k) - az(i-1,j,k)) * inv2dx[0];
                const Real dxay = (ay(i+1,j,k) - ay(i-1,j,k)) * inv2dx[0];
                cx(i,j,k) = -dzay;
                cy(i,j,k) = dzax - dxaz;
                cz(i,j,k) = dxay;
#else // WARPX_DIM_1D_Z
                // i = z; d/dx = d/dy = 0
                const Real dzay = (ay(i+1,j,k) - ay(i-1,j,k)) * inv2dx[2];
                const Real dzax = (ax(i+1,j,k) - ax(i-1,j,k)) * inv2dx[2];
                cx(i,j,k) = -dzay;
                cy(i,j,k) = dzax;
                cz(i,j,k) = 0.0_rt;
#endif
            });
    }
}

template <class T, class Ops>
void HybridPICPC<T,Ops>::Apply ( T& a_x, const T& a_b )
{
    using namespace amrex;

    BL_PROFILE("HybridPICPC::Apply()");

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(), "HybridPICPC::Apply() called on undefined object");

    auto& warpx = WarpX::GetInstance();
    const amrex::Periodicity& period = warpx.Geom(0).periodicity();

    for (int n = 0; n < 3; n++) {
        m_dE[n]->setVal(0.0);
        m_dW[n]->setVal(0.0);
    }

    const auto invdx2 = m_invdx2;
    const RT l2 = m_l2;
    const RT sigma = m_sigma;

    for (int sweep = 0; sweep < m_nsweeps; sweep++) {

        // --- Residuals of the split system (exact operators) ---
        if (sweep == 0) {
            // dE = dW = 0: r1 = b, r2 = 0
            for (int n = 0; n < 3; n++) {
                amrex::MultiFab::Copy(*m_r1[n], *a_b.getArrayVec()[0][n], 0, 0, 1, 0);
                m_r2[n]->setVal(0.0);
            }
        } else {
            for (int n = 0; n < 3; n++) {
                m_dE[n]->FillBoundary(period);
                m_dW[n]->FillBoundary(period);
            }
            // tmp = curl dE ; r2 = curl tmp - dW
            ComputeCurl(m_dE, m_tmp);
            for (int n = 0; n < 3; n++) { m_tmp[n]->FillBoundary(period); }
            ComputeCurl(m_tmp, m_r1); // borrow r1 as curlcurl scratch
            for (int n = 0; n < 3; n++) {
                amrex::MultiFab::Subtract(*m_r1[n], *m_dW[n], 0, 0, 1, 0);
                amrex::MultiFab::Copy(*m_r2[n], *m_r1[n], 0, 0, 1, 0);
            }
            // r1 = b - dE - [ cw (dW x bhat) + a_eta dW - a_H del2 dW ]
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*m_coef, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                Array4<const Real> const& coef = m_coef->const_array(mfi);
                Array4<const Real> const& bx0 = a_b.getArrayVec()[0][0]->const_array(mfi);
                Array4<const Real> const& by0 = a_b.getArrayVec()[0][1]->const_array(mfi);
                Array4<const Real> const& bz0 = a_b.getArrayVec()[0][2]->const_array(mfi);
                Array4<const Real> const& dEx = m_dE[0]->const_array(mfi);
                Array4<const Real> const& dEy = m_dE[1]->const_array(mfi);
                Array4<const Real> const& dEz = m_dE[2]->const_array(mfi);
                Array4<const Real> const& dWx = m_dW[0]->const_array(mfi);
                Array4<const Real> const& dWy = m_dW[1]->const_array(mfi);
                Array4<const Real> const& dWz = m_dW[2]->const_array(mfi);
                Array4<Real> const& r1x = m_r1[0]->array(mfi);
                Array4<Real> const& r1y = m_r1[1]->array(mfi);
                Array4<Real> const& r1z = m_r1[2]->array(mfi);
                amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        const Real bhx = coef(i,j,k,0);
                        const Real bhy = coef(i,j,k,1);
                        const Real bhz = coef(i,j,k,2);
                        const Real cw = coef(i,j,k,3);
                        const Real ae = coef(i,j,k,4);
                        const Real ah = coef(i,j,k,5);

                        // dW x bhat
                        const Real wxb_x = dWy(i,j,k)*bhz - dWz(i,j,k)*bhy;
                        const Real wxb_y = dWz(i,j,k)*bhx - dWx(i,j,k)*bhz;
                        const Real wxb_z = dWx(i,j,k)*bhy - dWy(i,j,k)*bhx;

                        // compact del2 dW (unused directions have invdx2 = 0)
                        Real lap_x = 0.0_rt, lap_y = 0.0_rt, lap_z = 0.0_rt;
#if defined(WARPX_DIM_3D)
                        lap_x = (dWx(i+1,j,k) - 2.0_rt*dWx(i,j,k) + dWx(i-1,j,k))*invdx2[0]
                              + (dWx(i,j+1,k) - 2.0_rt*dWx(i,j,k) + dWx(i,j-1,k))*invdx2[1]
                              + (dWx(i,j,k+1) - 2.0_rt*dWx(i,j,k) + dWx(i,j,k-1))*invdx2[2];
                        lap_y = (dWy(i+1,j,k) - 2.0_rt*dWy(i,j,k) + dWy(i-1,j,k))*invdx2[0]
                              + (dWy(i,j+1,k) - 2.0_rt*dWy(i,j,k) + dWy(i,j-1,k))*invdx2[1]
                              + (dWy(i,j,k+1) - 2.0_rt*dWy(i,j,k) + dWy(i,j,k-1))*invdx2[2];
                        lap_z = (dWz(i+1,j,k) - 2.0_rt*dWz(i,j,k) + dWz(i-1,j,k))*invdx2[0]
                              + (dWz(i,j+1,k) - 2.0_rt*dWz(i,j,k) + dWz(i,j-1,k))*invdx2[1]
                              + (dWz(i,j,k+1) - 2.0_rt*dWz(i,j,k) + dWz(i,j,k-1))*invdx2[2];
#elif defined(WARPX_DIM_XZ)
                        lap_x = (dWx(i+1,j,k) - 2.0_rt*dWx(i,j,k) + dWx(i-1,j,k))*invdx2[0]
                              + (dWx(i,j+1,k) - 2.0_rt*dWx(i,j,k) + dWx(i,j-1,k))*invdx2[2];
                        lap_y = (dWy(i+1,j,k) - 2.0_rt*dWy(i,j,k) + dWy(i-1,j,k))*invdx2[0]
                              + (dWy(i,j+1,k) - 2.0_rt*dWy(i,j,k) + dWy(i,j-1,k))*invdx2[2];
                        lap_z = (dWz(i+1,j,k) - 2.0_rt*dWz(i,j,k) + dWz(i-1,j,k))*invdx2[0]
                              + (dWz(i,j+1,k) - 2.0_rt*dWz(i,j,k) + dWz(i,j-1,k))*invdx2[2];
#else
                        lap_x = (dWx(i+1,j,k) - 2.0_rt*dWx(i,j,k) + dWx(i-1,j,k))*invdx2[2];
                        lap_y = (dWy(i+1,j,k) - 2.0_rt*dWy(i,j,k) + dWy(i-1,j,k))*invdx2[2];
                        lap_z = (dWz(i+1,j,k) - 2.0_rt*dWz(i,j,k) + dWz(i-1,j,k))*invdx2[2];
#endif
                        r1x(i,j,k) = bx0(i,j,k) - dEx(i,j,k)
                                   - (cw*wxb_x + ae*dWx(i,j,k) - ah*lap_x);
                        r1y(i,j,k) = by0(i,j,k) - dEy(i,j,k)
                                   - (cw*wxb_y + ae*dWy(i,j,k) - ah*lap_y);
                        r1z(i,j,k) = bz0(i,j,k) - dEz(i,j,k)
                                   - (cw*wxb_z + ae*dWz(i,j,k) - ah*lap_z);
                    });
            }
        }

        // --- Damped per-node block solve:
        //     dE += sigma * P^-1 (r1 - U r2),  dW += sigma * P^-1 (l2 r1 + r2) ---
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_coef, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<const Real> const& coef = m_coef->const_array(mfi);
            Array4<const Real> const& r1x = m_r1[0]->const_array(mfi);
            Array4<const Real> const& r1y = m_r1[1]->const_array(mfi);
            Array4<const Real> const& r1z = m_r1[2]->const_array(mfi);
            Array4<const Real> const& r2x = m_r2[0]->const_array(mfi);
            Array4<const Real> const& r2y = m_r2[1]->const_array(mfi);
            Array4<const Real> const& r2z = m_r2[2]->const_array(mfi);
            Array4<Real> const& dEx = m_dE[0]->array(mfi);
            Array4<Real> const& dEy = m_dE[1]->array(mfi);
            Array4<Real> const& dEz = m_dE[2]->array(mfi);
            Array4<Real> const& dWx = m_dW[0]->array(mfi);
            Array4<Real> const& dWy = m_dW[1]->array(mfi);
            Array4<Real> const& dWz = m_dW[2]->array(mfi);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    const Real bhx = coef(i,j,k,0);
                    const Real bhy = coef(i,j,k,1);
                    const Real bhz = coef(i,j,k,2);
                    const Real cw = coef(i,j,k,3);
                    const Real aU = coef(i,j,k,4) + coef(i,j,k,5)*l2;

                    const Real alpha = 1.0_rt + l2*aU;
                    const Real beta = l2*cw;
                    const Real inv_a2b2 = 1.0_rt / (alpha*alpha + beta*beta);

                    // s1 = r1 - U r2, U v = aU v + cw (v x bhat)
                    const Real r2xb_x = r2y(i,j,k)*bhz - r2z(i,j,k)*bhy;
                    const Real r2xb_y = r2z(i,j,k)*bhx - r2x(i,j,k)*bhz;
                    const Real r2xb_z = r2x(i,j,k)*bhy - r2y(i,j,k)*bhx;
                    const Real s1x = r1x(i,j,k) - (aU*r2x(i,j,k) + cw*r2xb_x);
                    const Real s1y = r1y(i,j,k) - (aU*r2y(i,j,k) + cw*r2xb_y);
                    const Real s1z = r1z(i,j,k) - (aU*r2z(i,j,k) + cw*r2xb_z);

                    // s2 = l2 r1 + r2
                    const Real s2x = l2*r1x(i,j,k) + r2x(i,j,k);
                    const Real s2y = l2*r1y(i,j,k) + r2y(i,j,k);
                    const Real s2z = l2*r1z(i,j,k) + r2z(i,j,k);

                    // P^-1 v = [alpha v - beta (v x bhat) + (beta^2/alpha)(v.bhat) bhat]/(alpha^2+beta^2)
                    const Real b2oa = beta*beta/alpha;

                    const Real s1xb_x = s1y*bhz - s1z*bhy;
                    const Real s1xb_y = s1z*bhx - s1x*bhz;
                    const Real s1xb_z = s1x*bhy - s1y*bhx;
                    const Real s1db = s1x*bhx + s1y*bhy + s1z*bhz;
                    dEx(i,j,k) += sigma*inv_a2b2*(alpha*s1x - beta*s1xb_x + b2oa*s1db*bhx);
                    dEy(i,j,k) += sigma*inv_a2b2*(alpha*s1y - beta*s1xb_y + b2oa*s1db*bhy);
                    dEz(i,j,k) += sigma*inv_a2b2*(alpha*s1z - beta*s1xb_z + b2oa*s1db*bhz);

                    const Real s2xb_x = s2y*bhz - s2z*bhy;
                    const Real s2xb_y = s2z*bhx - s2x*bhz;
                    const Real s2xb_z = s2x*bhy - s2y*bhx;
                    const Real s2db = s2x*bhx + s2y*bhy + s2z*bhz;
                    dWx(i,j,k) += sigma*inv_a2b2*(alpha*s2x - beta*s2xb_x + b2oa*s2db*bhx);
                    dWy(i,j,k) += sigma*inv_a2b2*(alpha*s2y - beta*s2xb_y + b2oa*s2db*bhy);
                    dWz(i,j,k) += sigma*inv_a2b2*(alpha*s2z - beta*s2xb_z + b2oa*s2db*bhz);
                });
        }
    }

    // x = dE
    for (int n = 0; n < 3; n++) {
        amrex::MultiFab::Copy(*a_x.getArrayVec()[0][n], *m_dE[n], 0, 0, 1, 0);
    }
}

#endif
