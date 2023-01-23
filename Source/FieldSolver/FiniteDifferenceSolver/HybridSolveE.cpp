#include "FiniteDifferenceSolver.H"

#ifdef WARPX_DIM_RZ
#   include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
#   include "FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#endif
#include "HybridModel/HybridModel.H"
#include "Utils/CoarsenIO.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

using namespace amrex;

void FiniteDifferenceSolver::HybridSolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield_old,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, std::unique_ptr<HybridModel> const& hybrid_model,
    DtType a_dt_type )
{

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
    // WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
    //     !m_do_nodal, "hybrid E-solve does not work for nodal");

    if (m_fdtd_algo == ElectromagneticSolverAlgo::Hybrid) {
#ifdef WARPX_DIM_RZ
        HybridSolveECylindrical <CylindricalYeeAlgorithm> (
            Efield, Jfield, Bfield, Jifield, Jifield_old, rhofield,
            edge_lengths, lev, hybrid_model, a_dt_type
        );
#else
        if (m_do_nodal) {
            HybridSolveECartesian <CartesianNodalAlgorithm> (
                Efield, Jfield, Bfield, Jifield, Jifield_old, rhofield,
                edge_lengths, lev, hybrid_model, a_dt_type
            );
        }
        else {
            HybridSolveECartesian <CartesianYeeAlgorithm> (
                Efield, Jfield, Bfield, Jifield, Jifield_old, rhofield,
                edge_lengths, lev, hybrid_model, a_dt_type
            );
        }
#endif
    } else {
        amrex::Abort(Utils::TextMsg::Err(
            "HybridSolveE: The hybrid electromagnetic solver algorithm must be used"));
    }
}

#ifdef WARPX_DIM_RZ
template<typename T_Algo>
void FiniteDifferenceSolver::HybridSolveECylindrical (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield_old,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, std::unique_ptr<HybridModel> const& hybrid_model,
    DtType a_dt_type )
{
#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif
    amrex::ignore_unused(Efield, Bfield, Jifield, rhofield, edge_lengths);
    amrex::Abort(Utils::TextMsg::Err(
        "currently hybrid E-solve does not work for RZ"));
}

#else

template<typename T_Algo>
void FiniteDifferenceSolver::HybridSolveECartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield_old,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, std::unique_ptr<HybridModel> const& hybrid_model,
    DtType a_dt_type )
{
#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif

    // for the profiler
    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // get hybrid model parameters
    auto n0 = hybrid_model->m_n0_ref;
    auto T0 = hybrid_model->m_elec_temp;
    auto gamma = hybrid_model->m_gamma;
    auto eta = hybrid_model->m_eta;

    // Index type required for calling CoarsenIO::Interp to interpolate fields
    // from their respective staggering to the Ex, Ey, Ez locations
    amrex::GpuArray<int, 3> const& rho_stag = hybrid_model->rho_IndexType;
    amrex::GpuArray<int, 3> const& Jx_stag = hybrid_model->Jx_IndexType;
    amrex::GpuArray<int, 3> const& Jy_stag = hybrid_model->Jy_IndexType;
    amrex::GpuArray<int, 3> const& Jz_stag = hybrid_model->Jz_IndexType;
    amrex::GpuArray<int, 3> const& Bx_stag = hybrid_model->Bx_IndexType;
    amrex::GpuArray<int, 3> const& By_stag = hybrid_model->By_IndexType;
    amrex::GpuArray<int, 3> const& Bz_stag = hybrid_model->Bz_IndexType;
    amrex::GpuArray<int, 3> const& Ex_stag = hybrid_model->Ex_IndexType;
    amrex::GpuArray<int, 3> const& Ey_stag = hybrid_model->Ey_IndexType;
    amrex::GpuArray<int, 3> const& Ez_stag = hybrid_model->Ez_IndexType;

    // The "coarsening is just 1 i.e. no coarsening"
    amrex::GpuArray<int, 3> const& coarsen = {1, 1, 1};

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        // Extract field data for this grid/tile
        Array4<Real> const& Ex = Efield[0]->array(mfi);
        Array4<Real> const& Ey = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);
        Array4<Real> const& Jx = Jfield[0]->array(mfi);
        Array4<Real> const& Jy = Jfield[1]->array(mfi);
        Array4<Real> const& Jz = Jfield[2]->array(mfi);
        Array4<Real> const& Bx = Bfield[0]->array(mfi);
        Array4<Real> const& By = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);
        Array4<Real> const& Jix = Jifield[0]->array(mfi);
        Array4<Real> const& Jiy = Jifield[1]->array(mfi);
        Array4<Real> const& Jiz = Jifield[2]->array(mfi);
        Array4<Real> const& rho = rhofield->array(mfi);

        Array4<Real> const& Jix_old = Jifield_old[0]->array(mfi);
        Array4<Real> const& Jiy_old = Jifield_old[1]->array(mfi);
        Array4<Real> const& Jiz_old = Jifield_old[2]->array(mfi);

#ifdef AMREX_USE_EB
        amrex::Array4<amrex::Real> const& lx = edge_lengths[0]->array(mfi);
        amrex::Array4<amrex::Real> const& ly = edge_lengths[1]->array(mfi);
        amrex::Array4<amrex::Real> const& lz = edge_lengths[2]->array(mfi);
#endif

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        int const n_coefs_x = m_stencil_coefs_x.size();
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        int const n_coefs_y = m_stencil_coefs_y.size();
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        int const n_coefs_z = m_stencil_coefs_z.size();

        // Extract tileboxes for which to loop
        Box const& tjx  = mfi.tilebox(Jfield[0]->ixType().toIntVect());
        Box const& tjy  = mfi.tilebox(Jfield[1]->ixType().toIntVect());
        Box const& tjz  = mfi.tilebox(Jfield[2]->ixType().toIntVect());

        // First calculate the total current using Ampere's law on the
        // same grid as the E-field
        amrex::ParallelFor(tjx, tjy, tjz,

            // Jx calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip if this cell is fully covered by embedded boundaries
                if (lx(i, j, k) <= 0) return;
#endif
                Jx(i, j, k) = (
                    - T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k)
                    + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k)
                ) / PhysConst::mu0;
            },

            // Jy calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip field push if this cell is fully covered by embedded boundaries
#ifdef WARPX_DIM_3D
                if (ly(i,j,k) <= 0) return;
#elif defined(WARPX_DIM_XZ)
                // In XZ Jy is associated with a mesh node, so we need to check if the mesh node is covered
                amrex::ignore_unused(ly);
                if (lx(i, j, k)<=0 || lx(i-1, j, k)<=0 || lz(i, j-1, k)<=0 || lz(i, j, k)<=0) return;
#endif
#endif
                Jy(i, j, k) = (
                    - T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k)
                    + T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k)
                ) / PhysConst::mu0;
            },

            // Jz calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip if this cell is fully covered by embedded boundaries
                if (lz(i,j,k) <= 0) return;
#endif
                Jz(i, j, k) = (
                    - T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k)
                    + T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k)
                ) / PhysConst::mu0;
            }
        );

        Box const& tex  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tey  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tez  = mfi.tilebox(Efield[2]->ixType().toIntVect());

        // Loop over the cells and update the E field
        amrex::ParallelFor(tex, tey, tez,

            // Ex calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip if this cell is fully covered by embedded boundaries
                if (lx(i, j, k) <= 0) return;
#endif

                // allocate variables for all field quantities (interpolated onto the correct grid)
                Real rho_interp, jey_interp, jez_interp, jy_interp, jz_interp, grad_p;

                // get the appropriate charge density
                if (a_dt_type == DtType::Full) {
                    // use rho^{n}
                    rho_interp = CoarsenIO::Interp(rho, rho_stag, Ex_stag, coarsen, i, j, k, 0);
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use rho^{n+1/2}
                    rho_interp = 0.5_rt * (
                        CoarsenIO::Interp(rho, rho_stag, Ex_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(rho, rho_stag, Ex_stag, coarsen, i, j, k, 1)
                    );
                } else if (a_dt_type == DtType::SecondHalf) {
                    // use rho^{n+1}
                    rho_interp = CoarsenIO::Interp(rho, rho_stag, Ex_stag, coarsen, i, j, k, 1);
                }

                if (rho_interp == 0._rt) {
                    Ex(i, j, k) = 0._rt;
                    return;
                }
                // rho_interp = n0 * PhysConst::q_e;

                // interpolate the total current to the appropriate grid
                jy_interp = CoarsenIO::Interp(Jy, Jy_stag, Ex_stag, coarsen, i, j, k, 0);
                jz_interp = CoarsenIO::Interp(Jz, Jz_stag, Ex_stag, coarsen, i, j, k, 0);

                // get the electron current at the appropriate time and location
                if (a_dt_type == DtType::Full) {
                    // use J^{n}
                    jey_interp = jy_interp - 0.5_rt * (
                        CoarsenIO::Interp(Jiy_old, Jy_stag, Ex_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(Jiy, Jy_stag, Ex_stag, coarsen, i, j, k, 0)
                    );
                    jez_interp = jz_interp - 0.5_rt * (
                        CoarsenIO::Interp(Jiz_old, Jz_stag, Ex_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(Jiz, Jz_stag, Ex_stag, coarsen, i, j, k, 0)
                    );
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use J^{n+1/2}
                    jey_interp = jy_interp - CoarsenIO::Interp(Jiy, Jy_stag, Ex_stag, coarsen, i, j, k, 0);
                    jez_interp = jz_interp - CoarsenIO::Interp(Jiz, Jz_stag, Ex_stag, coarsen, i, j, k, 0);
                }
                else if (a_dt_type == DtType::SecondHalf) {
                    // use J^{n+1}
                    jey_interp = jy_interp - 0.5_rt * (
                        3._rt * CoarsenIO::Interp(Jiy, Jy_stag, Ex_stag, coarsen, i, j, k, 0)
                        - CoarsenIO::Interp(Jiy_old, Jy_stag, Ex_stag, coarsen, i, j, k, 0)
                    );
                    jez_interp = jz_interp - 0.5_rt * (
                        3._rt * CoarsenIO::Interp(Jiz, Jz_stag, Ex_stag, coarsen, i, j, k, 0)
                        - CoarsenIO::Interp(Jiz_old, Jz_stag, Ex_stag, coarsen, i, j, k, 0)
                    );
                }

                auto const By_interp = CoarsenIO::Interp(
                    By, By_stag, Ex_stag, coarsen, i, j, k, 0
                );
                auto const Bz_interp = CoarsenIO::Interp(
                    Bz, Bz_stag, Ex_stag, coarsen, i, j, k, 0
                );

#if (defined WARPX_DIM_1D_Z)
                grad_p = 0._rt; // 1D Cartesian: derivative along x is 0
#else
                if (a_dt_type == DtType::Full) {
                    // use rho^{n}
                    grad_p = coefs_x[0]*(
                        ElectronPressure::get_pressure(n0, T0, gamma, rho(i+1, j, k, 0))
                        - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 0))
                    );
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use rho^{n+1/2}
                    grad_p = coefs_x[0]*(
                        ElectronPressure::get_pressure(
                            n0, T0, gamma,
                            0.5_rt * (rho(i+1, j, k, 1) + rho(i+1, j, k, 0))
                        )
                        - ElectronPressure::get_pressure(
                            n0, T0, gamma,
                            0.5_rt * (rho(i, j, k, 1) + rho(i, j, k, 0))
                        )
                    );
                } else if (a_dt_type == DtType::SecondHalf) {
                    // use rho^{n+1}
                    grad_p = coefs_x[0]*(
                        ElectronPressure::get_pressure(n0, T0, gamma, rho(i+1, j, k, 1))
                        - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 1))
                    );
                }
#endif

                // calculate the full current
                // auto const Jx =  T_Algo::Jx(
                //         By, Bz, Bx_stag, By_stag, coefs_y, coefs_z, Ex_stag, i, j, k
                // );
                // auto const Jy =  T_Algo::Jy(
                //         Bx, Bz, Bx_stag, Bz_stag, coefs_x, coefs_z, Ex_stag, i, j, k
                // );
                // auto const Jz = T_Algo::Jz(
                //         Bx, By, Bx_stag, By_stag, coefs_x, coefs_y, Ex_stag, i, j, k
                // );
                // auto const Jx = 0._rt;
                // auto const Jy = 0._rt;
                // auto const Jz = 0._rt;

                Ex(i, j, k) = (
                    jey_interp * Bz_interp - jez_interp * By_interp
                    - grad_p
                ) / rho_interp + eta * Jx(i, j, k);

                // if (i < 100) {
                //     amrex::Print() << Ex(i,j,k) << "  "
                //         << Jx << "   " <<  Jy << "  " <<  Jz << "  "
                //         << jz_interp << "   " <<  jy_interp << "  "
                //         << Bz_interp
                //     << std::endl;
                // }
            },

            // Ey calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip field solve if this cell is fully covered by embedded boundaries
#ifdef WARPX_DIM_3D
                if (ly(i,j,k) <= 0) return;
#elif defined(WARPX_DIM_XZ)
                //In XZ Ey is associated with a mesh node, so we need to check if the mesh node is covered
                amrex::ignore_unused(ly);
                if (lx(i, j, k)<=0 || lx(i-1, j, k)<=0 || lz(i, j-1, k)<=0 || lz(i, j, k)<=0) return;
#endif
#endif

                // allocate variables for all interpolated (onto correct grid) field quantities
                Real rho_interp, jex_interp, jez_interp, jx_interp, jz_interp, grad_p;

                // get the appropriate charge density
                if (a_dt_type == DtType::Full) {
                    // use rho^{n}
                    rho_interp = CoarsenIO::Interp(rho, rho_stag, Ey_stag, coarsen, i, j, k, 0);
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use rho^{n+1/2}
                    rho_interp = 0.5_rt * (
                        CoarsenIO::Interp(rho, rho_stag, Ey_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(rho, rho_stag, Ey_stag, coarsen, i, j, k, 1)
                    );
                }
                else if (a_dt_type == DtType::SecondHalf) {
                    // use rho^{n+1}
                    rho_interp = CoarsenIO::Interp(rho, rho_stag, Ey_stag, coarsen, i, j, k, 1);
                }

                if (rho_interp == 0._rt) {
                    Ey(i, j, k) = 0._rt;
                    return;
                }
                // rho_interp = n0 * PhysConst::q_e;

                // interpolate the total current to the appropriate grid
                jx_interp = CoarsenIO::Interp(Jx, Jx_stag, Ey_stag, coarsen, i, j, k, 0);
                jz_interp = CoarsenIO::Interp(Jz, Jz_stag, Ey_stag, coarsen, i, j, k, 0);

                // get the electron current at the appropriate time and location
                if (a_dt_type == DtType::Full) {
                    // use J^{n}
                    jex_interp = jx_interp - 0.5_rt * (
                        CoarsenIO::Interp(Jix_old, Jx_stag, Ey_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(Jix, Jx_stag, Ey_stag, coarsen, i, j, k, 0)
                    );
                    jez_interp = jz_interp - 0.5_rt * (
                        CoarsenIO::Interp(Jiz_old, Jz_stag, Ey_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(Jiz, Jz_stag, Ey_stag, coarsen, i, j, k, 0)
                    );
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use J^{n+1/2}
                    jex_interp = jx_interp - CoarsenIO::Interp(Jix, Jx_stag, Ey_stag, coarsen, i, j, k, 0);
                    jez_interp = jz_interp - CoarsenIO::Interp(Jiz, Jz_stag, Ey_stag, coarsen, i, j, k, 0);
                } else if (a_dt_type == DtType::SecondHalf) {
                    // use J^{n+1}
                    jex_interp = jx_interp - 0.5_rt * (
                        3._rt * CoarsenIO::Interp(Jix, Jx_stag, Ey_stag, coarsen, i, j, k, 0)
                        - CoarsenIO::Interp(Jix_old, Jx_stag, Ey_stag, coarsen, i, j, k, 0)
                    );
                    jez_interp = jz_interp - 0.5_rt * (
                        3._rt * CoarsenIO::Interp(Jiz, Jz_stag, Ey_stag, coarsen, i, j, k, 0)
                        - CoarsenIO::Interp(Jiz_old, Jz_stag, Ey_stag, coarsen, i, j, k, 0)
                    );
                }

                auto const Bx_interp = CoarsenIO::Interp(
                    Bx, Bx_stag, Ey_stag, coarsen, i, j, k, 0
                );
                auto const Bz_interp = CoarsenIO::Interp(
                    Bz, Bz_stag, Ey_stag, coarsen, i, j, k, 0
                );

#if (defined WARPX_DIM_1D_Z) || (defined WARPX_DIM_XZ)
                grad_p = 0._rt; // 1D and 2D Cartesian: derivative along y is 0
#else
                if (a_dt_type == DtType::Full) {
                    // use rho^{n}
                    grad_p = coefs_y[0]*(
                        ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j+1, k, 0))
                        - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 0))
                    );
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use rho^{n+1/2}
                    grad_p = coefs_y[0]*(
                        ElectronPressure::get_pressure(
                            n0, T0, gamma,
                            0.5_rt * (rho(i, j+1, k, 1) + rho(i, j+1, k, 0))
                        )
                        - ElectronPressure::get_pressure(
                            n0, T0, gamma,
                            0.5_rt * (rho(i, j, k, 1) + rho(i, j, k, 0))
                        )
                    );
                } else if (a_dt_type == DtType::SecondHalf) {
                    // use rho^{n+1}
                    grad_p = coefs_y[0]*(
                        ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j+1, k, 1))
                        - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 1))
                    );
                }
#endif

                // calculate the full current
                // auto const Jx =  T_Algo::Jx(
                //         By, Bz, By_stag, Bz_stag, coefs_y, coefs_z, Ey_stag, i, j, k
                // );
                // auto const Jy =  T_Algo::Jy(
                //         Bx, Bz, Bx_stag, Bz_stag, coefs_x, coefs_z, Ey_stag, i, j, k
                // );
                // auto const Jz = T_Algo::Jz(
                //         Bx, By, Bx_stag, By_stag, coefs_x, coefs_y, Ey_stag, i, j, k
                // );
                // auto const Jx = 0._rt;
                // auto const Jy = 0._rt;
                // auto const Jz = 0._rt;

                Ey(i, j, k) = (
                    jez_interp * Bx_interp - jex_interp * Bz_interp
                    - grad_p
                ) / rho_interp + eta * Jy(i, j, k);
            },

            // Ez calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){

#ifdef AMREX_USE_EB
                // Skip field solve if this cell is fully covered by embedded boundaries
                if (lz(i,j,k) <= 0) return;
#endif

                // allocate variables for all interpolated (onto correct grid) field quantities
                Real rho_interp, jex_interp, jey_interp, jx_interp, jy_interp, grad_p;

                // get the appropriate charge density
                if (a_dt_type == DtType::Full) {
                    // use rho^{n}
                    rho_interp = CoarsenIO::Interp(rho, rho_stag, Ez_stag, coarsen, i, j, k, 0);
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use rho^{n+1/2}
                    rho_interp = 0.5_rt * (
                        CoarsenIO::Interp(rho, rho_stag, Ez_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(rho, rho_stag, Ez_stag, coarsen, i, j, k, 1)
                    );
                } else if (a_dt_type == DtType::SecondHalf) {
                    // use rho^{n+1}
                    rho_interp = CoarsenIO::Interp(rho, rho_stag, Ez_stag, coarsen, i, j, k, 1);
                }

                if (rho_interp == 0._rt) {
                    Ez(i, j, k) = 0._rt;
                    return;
                }
                // rho_interp = n0 * PhysConst::q_e;

                // interpolate the total current to the appropriate grid
                jx_interp = CoarsenIO::Interp(Jx, Jx_stag, Ez_stag, coarsen, i, j, k, 0);
                jy_interp = CoarsenIO::Interp(Jy, Jy_stag, Ez_stag, coarsen, i, j, k, 0);

                // get the electron current at the appropriate time and location
                if (a_dt_type == DtType::Full) {
                    // use J^{n}
                    jex_interp = jx_interp - 0.5_rt * (
                        CoarsenIO::Interp(Jix_old, Jx_stag, Ez_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(Jix, Jx_stag, Ez_stag, coarsen, i, j, k, 0)
                    );
                    jey_interp = jy_interp - 0.5_rt * (
                        CoarsenIO::Interp(Jiy_old, Jy_stag, Ez_stag, coarsen, i, j, k, 0)
                        + CoarsenIO::Interp(Jiy, Jy_stag, Ez_stag, coarsen, i, j, k, 0)
                    );
                } else if (a_dt_type == DtType::FirstHalf) {
                    // use J^{n+1/2}
                    jex_interp = jx_interp - CoarsenIO::Interp(Jix, Jx_stag, Ez_stag, coarsen, i, j, k, 0);
                    jey_interp = jy_interp - CoarsenIO::Interp(Jiy, Jy_stag, Ez_stag, coarsen, i, j, k, 0);
                }
                else if (a_dt_type == DtType::SecondHalf) {
                    // use J^{n+1}
                    jex_interp = jx_interp - 0.5_rt * (
                        3._rt * CoarsenIO::Interp(Jix, Jx_stag, Ez_stag, coarsen, i, j, k, 0)
                        - CoarsenIO::Interp(Jix_old, Jx_stag, Ez_stag, coarsen, i, j, k, 0)
                    );
                    jey_interp = jy_interp - 0.5_rt * (
                        3._rt * CoarsenIO::Interp(Jiy, Jy_stag, Ez_stag, coarsen, i, j, k, 0)
                        - CoarsenIO::Interp(Jiy_old, Jy_stag, Ez_stag, coarsen, i, j, k, 0)
                    );
                }

                auto const Bx_interp = CoarsenIO::Interp(
                    Bx, Bx_stag, Ez_stag, coarsen, i, j, k, 0
                );
                auto const By_interp = CoarsenIO::Interp(
                    By, By_stag, Ez_stag, coarsen, i, j, k, 0
                );

#if (defined WARPX_DIM_1D_Z)
                auto i1 = i+1;
                auto j1 = j;
                auto k1 = k;
#elif (defined WARPX_DIM_XZ)
                auto i1 = i;
                auto j1 = j+1;
                auto k1 = k;
#else
                auto i1 = i;
                auto j1 = j;
                auto k1 = k+1;
#endif
                if (a_dt_type == DtType::Full) {
                    // use rho^{n}
                    grad_p = coefs_z[0]*(
                        ElectronPressure::get_pressure(n0, T0, gamma, rho(i1, j1, k1, 0))
                        - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 0))
                    );
                }
                else if (a_dt_type == DtType::FirstHalf) {
                    // use rho^{n+1/2}
                    grad_p = coefs_z[0]*(
                        ElectronPressure::get_pressure(
                            n0, T0, gamma,
                            0.5_rt * (rho(i1, j1, k1, 1) + rho(i1, j1, k1, 0))
                        )
                        - ElectronPressure::get_pressure(
                            n0, T0, gamma,
                            0.5_rt * (rho(i, j, k, 1) + rho(i, j, k, 0))
                        )
                    );
                }
                else if (a_dt_type == DtType::SecondHalf) {
                    // use rho^{n+1}
                    grad_p = coefs_z[0]*(
                        ElectronPressure::get_pressure(n0, T0, gamma, rho(i1, j1, k1, 1))
                        - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 1))
                    );
                }

                // calculate the full current
                // auto const Jx =  T_Algo::Jx(
                //         By, Bz, By_stag, Bz_stag, coefs_y, coefs_z, Ez_stag, i, j, k
                // );
                // auto const Jy = T_Algo::Jy(
                //         Bx, Bz, Bx_stag, Bz_stag, coefs_x, coefs_z, Ez_stag, i, j, k
                // );
                // auto const Jz = T_Algo::Jz(
                //         Bx, By, Bx_stag, By_stag, coefs_x, coefs_y, Ez_stag, i, j, k
                // );
                // auto const Jx = 0._rt;
                // auto const Jy = 0._rt;
                // auto const Jz = 0._rt;

                Ez(i, j, k) = (
                    jex_interp * By_interp - jey_interp * Bx_interp
                    - grad_p
                ) / rho_interp + eta * Jz(i, j, k);

                // amrex::Print() << "[ " << i << ", " << j << "]  " << Ez(i, j, k) << "  " << rho_interp << "  " << Jy << "  " << grad_p << "  " << Bx_interp << std::endl;

            }

        );

        // // If F is not a null pointer, further update E using the grad(F) term
        // // (hyperbolic correction for errors in charge conservation)
        // if (Ffield) {

        //     // Extract field data for this grid/tile
        //     Array4<Real> F = Ffield->array(mfi);

        //     // Loop over the cells and update the fields
        //     amrex::ParallelFor(tex, tey, tez,

        //         [=] AMREX_GPU_DEVICE (int i, int j, int k){
        //             Ex(i, j, k) += c2 * dt * T_Algo::UpwardDx(F, coefs_x, n_coefs_x, i, j, k);
        //         },
        //         [=] AMREX_GPU_DEVICE (int i, int j, int k){
        //             Ey(i, j, k) += c2 * dt * T_Algo::UpwardDy(F, coefs_y, n_coefs_y, i, j, k);
        //         },
        //         [=] AMREX_GPU_DEVICE (int i, int j, int k){
        //             Ez(i, j, k) += c2 * dt * T_Algo::UpwardDz(F, coefs_z, n_coefs_z, i, j, k);
        //         }

        //     );

        // }

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }

}
#endif