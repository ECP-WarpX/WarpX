#include "Utils/WarpXAlgorithmSelection.H"
#include "FiniteDifferenceSolver.H"
#ifdef WARPX_DIM_RZ
    // currently works only for 3D
#else
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#endif
#include "Utils/WarpXConst.H"
#include <AMReX_Gpu.H>
#include <WarpX.H>

using namespace amrex;

void FiniteDifferenceSolver::MacroscopicEvolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    amrex::Real const dt, std::unique_ptr<MacroscopicProperties> const& macroscopic_properties ) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    amrex::Abort("currently macro E-push does not work for RZ");
#else
    if (m_do_nodal) {
        amrex::Abort(" macro E-push does not work for nodal ");

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {

        MacroscopicEvolveECartesian <CartesianYeeAlgorithm> ( Efield, Bfield, Jfield, dt,
                                                              macroscopic_properties );

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        // Note : EvolveE is the same for CKC and Yee.
        // In the templated Yee and CKC calls, the core operations for EvolveE is the same.
        MacroscopicEvolveECartesian <CartesianCKCAlgorithm> ( Efield, Bfield, Jfield, dt,
                                                              macroscopic_properties );
    } else {
        amrex::Abort("Unknown algorithm");
    }
#endif

}


#ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::MacroscopicEvolveECartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    amrex::Real const dt, std::unique_ptr<MacroscopicProperties> const& macroscopic_properties ) {

    const int &macroscopic_solver_algo = WarpX::macroscopic_solver_algo;
    Real sigma = macroscopic_properties->sigma();
    Real const mu = macroscopic_properties->mu();
    Real const epsilon = macroscopic_properties->epsilon();

    Real alpha = 0._rt;
    Real beta  = 0._rt;
    Real fac1 = 0._rt;
    Real inv_fac = 0._rt;

    if (macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler) {
        fac1 = sigma * dt / epsilon;
        inv_fac = 1._rt / ( 1._rt + fac1);
        alpha = inv_fac;
        beta  = dt * inv_fac / epsilon;
    } else if (macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff) {
        fac1 = 0.5_rt * sigma * dt / epsilon;
        inv_fac = 1._rt / ( 1._rt + fac1);
        alpha = (1.0_rt - fac1) * inv_fac;
        beta  = dt * inv_fac / epsilon;
    }

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        Array4<Real> const& Ex = Efield[0]->array(mfi);
        Array4<Real> const& Ey = Efield[1]->array(mfi);
        Array4<Real> const& Ez = Efield[2]->array(mfi);
        Array4<Real> const& Bx = Bfield[0]->array(mfi);
        Array4<Real> const& By = Bfield[1]->array(mfi);
        Array4<Real> const& Bz = Bfield[2]->array(mfi);

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
                Ex(i, j, k) = alpha * Ex(i, j, k) + (beta/mu)
                     * ( - T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k)
                         + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k));
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ey(i, j, k) = alpha * Ey(i, j, k) + (beta/mu)
                     * ( - T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k)
                         + T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k));
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Ez(i, j, k) = alpha * Ez(i, j, k) + (beta/mu)
                     * ( - T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k)
                         + T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k));
            }

        );

        // update E using J, if source currents are specified.
        if (Jfield[0]) {
            Array4<Real> const& jx = Jfield[0]->array(mfi);
            Array4<Real> const& jy = Jfield[1]->array(mfi);
            Array4<Real> const& jz = Jfield[2]->array(mfi);

            amrex::ParallelFor(tex, tey, tez,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ex(i, j, k) += -beta * jx(i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ey(i, j, k) += -beta * jy(i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Ez(i, j, k) += -beta * jz(i, j, k);
                }
            );
        }
    }

}

#endif // corresponds to ifndef WARPX_DIM_RZ
