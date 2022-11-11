#include "FiniteDifferenceSolver.H"

#ifdef WARPX_DIM_RZ
#   include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
#   include "HybridModel/CartesianHybridYeeAlgorithm.H"
#endif
#include "HybridModel/HybridModel.H"
#include "Utils/CoarsenIO.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

using namespace amrex;

void FiniteDifferenceSolver::HybridSolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3>& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, std::unique_ptr<HybridModel> const& hybrid_model )
{

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_do_nodal, "hybrid E-solve does not work for nodal");

    if (m_fdtd_algo == ElectromagneticSolverAlgo::Hybrid) {
#ifdef WARPX_DIM_RZ
        HybridSolveECylindrical <CylindricalYeeAlgorithm> (
            Efield, Bfield, Jfield, rhofield, edge_lengths, lev
        );
#else
        HybridSolveECartesian <CartesianHybridYeeAlgorithm> (
            Efield, Bfield, Jfield, rhofield, edge_lengths, lev, hybrid_model
        );
#endif
    } else {
        amrex::Abort(Utils::TextMsg::Err(
            "HybridSolveE: The hybrid electromagnetic solver algorithm must be used"));
    }
}

#ifdef WARPX_DIM_RZ
template<typename T_Algo>
void FiniteDifferenceSolver::HybridSolveECylindrical (
    std::array< std::unique_ptr<amrex::MultiFab>, 3>& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, std::unique_ptr<HybridModel> const& hybrid_model )
{
#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif
    amrex::ignore_unused(Efield, Bfield, Jfield, rhofield, edge_lengths);
    amrex::Abort(Utils::TextMsg::Err(
        "currently hybrid E-solve does not work for RZ"));
}

#else

template<typename T_Algo>
void FiniteDifferenceSolver::HybridSolveECartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3>& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, std::unique_ptr<HybridModel> const& hybrid_model )
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
        Array4<Real> const& Bx = Bfield[0]->array(mfi);
        Array4<Real> const& By = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);
        Array4<Real> const& jx = Jfield[0]->array(mfi);
        Array4<Real> const& jy = Jfield[1]->array(mfi);
        Array4<Real> const& jz = Jfield[2]->array(mfi);
        Array4<Real> const& rho = rhofield->array(mfi);

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
        Box const& tex  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tey  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tez  = mfi.tilebox(Efield[2]->ixType().toIntVect());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tex, tey, tez,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip field solve if this cell is fully covered by embedded boundaries
                if (lx(i, j, k) <= 0) return;
#endif

                auto const rho_interp = CoarsenIO::Interp(
                    rho, rho_stag, Ex_stag, coarsen, i, j, k, 0
                );
                auto const Jy_interp = CoarsenIO::Interp(
                    jy, Jy_stag, Ex_stag, coarsen, i, j, k, 0
                );
                auto const Jz_interp = CoarsenIO::Interp(
                    jz, Jz_stag, Ex_stag, coarsen, i, j, k, 0
                );
                auto const By_interp = CoarsenIO::Interp(
                    By, By_stag, Ex_stag, coarsen, i, j, k, 0
                );
                auto const Bz_interp = CoarsenIO::Interp(
                    Bz, Bz_stag, Ex_stag, coarsen, i, j, k, 0
                );

                // auto const curlB_y =  T_Algo::curlB_y_Ex(
                //         Bx, Bz, Bx_stag, Bz_stag, coefs_x, coefs_z, i, j, k
                // );
                // auto const curlB_z = T_Algo::curlB_z_Ex(
                //         Bx, By, Bx_stag, By_stag, coefs_x, coefs_y, i, j, k
                // );

                // TODO: should use average rho if middle step...
                auto const grad_p =
#if (defined WARPX_DIM_1D_Z)
                0._rt; // 1D Cartesian: derivative along x is 0
#else
                coefs_x[0]*(
                    ElectronPressure::get_pressure(n0, T0, gamma, rho(i+1, j, k, 1))
                    - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 1))
                );
#endif

                Ex(i, j, k) = (
                    - Bz_interp * (Jy_interp) // - curlB_y / PhysConst::mu0)
                    + By_interp * (Jz_interp) // - curlB_z / PhysConst::mu0)
                    - grad_p
                ) / rho_interp;
            },

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

                auto const rho_interp = CoarsenIO::Interp(
                    rho, rho_stag, Ey_stag, coarsen, i, j, k, 0
                );
                auto const Jx_interp = CoarsenIO::Interp(
                    jx, Jx_stag, Ey_stag, coarsen, i, j, k, 0
                );
                auto const Jz_interp = CoarsenIO::Interp(
                    jz, Jz_stag, Ey_stag, coarsen, i, j, k, 0
                );
                auto const Bx_interp = CoarsenIO::Interp(
                    Bx, Bx_stag, Ey_stag, coarsen, i, j, k, 0
                );
                auto const Bz_interp = CoarsenIO::Interp(
                    Bz, Bz_stag, Ey_stag, coarsen, i, j, k, 0
                );

                // auto const curlB_y =  T_Algo::curlB_y_Ex(
                //         Bx, Bz, Bx_stag, Bz_stag, coefs_x, coefs_z, i, j, k
                // );
                // auto const curlB_z = T_Algo::curlB_z_Ex(
                //         Bx, By, Bx_stag, By_stag, coefs_x, coefs_y, i, j, k
                // );

                // TODO: should use average rho if middle step...
                auto const grad_p =
#if (defined WARPX_DIM_1D_Z) || (defined WARPX_DIM_XZ)
                0._rt; // 1D Cartesian: derivative along y is 0
#else
                coefs_y[0]*(
                    ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j+1, k, 1))
                    - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 1))
                );
#endif

                Ey(i, j, k) = (
                    - Bz_interp * (Jx_interp) // - curlB_y / PhysConst::mu0)
                    + Bx_interp * (Jz_interp) // - curlB_z / PhysConst::mu0)
                    - grad_p
                ) / rho_interp;
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){

#ifdef AMREX_USE_EB
                // Skip field solve if this cell is fully covered by embedded boundaries
                if (lz(i,j,k) <= 0) return;
#endif

                auto const rho_interp = CoarsenIO::Interp(
                    rho, rho_stag, Ez_stag, coarsen, i, j, k, 0
                );
                auto const Jx_interp = CoarsenIO::Interp(
                    jx, Jx_stag, Ez_stag, coarsen, i, j, k, 0
                );
                auto const Jy_interp = CoarsenIO::Interp(
                    jy, Jy_stag, Ez_stag, coarsen, i, j, k, 0
                );
                auto const Bx_interp = CoarsenIO::Interp(
                    Bx, Bx_stag, Ez_stag, coarsen, i, j, k, 0
                );
                auto const By_interp = CoarsenIO::Interp(
                    By, By_stag, Ez_stag, coarsen, i, j, k, 0
                );

                // auto const curlB_y =  T_Algo::curlB_y_Ex(
                //         Bx, Bz, Bx_stag, Bz_stag, coefs_x, coefs_z, i, j, k
                // );
                // auto const curlB_z = T_Algo::curlB_z_Ex(
                //         Bx, By, Bx_stag, By_stag, coefs_x, coefs_y, i, j, k
                // );

                // TODO: should use average rho if middle step...
                auto const grad_p = coefs_z[0]*(
                    ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k+1, 1))
                    - ElectronPressure::get_pressure(n0, T0, gamma, rho(i, j, k, 1))
                );

                Ez(i, j, k) = (
                    - Bx_interp * (Jy_interp) // - curlB_y / PhysConst::mu0)
                    + By_interp * (Jx_interp) // - curlB_z / PhysConst::mu0)
                    - grad_p
                ) / rho_interp;
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