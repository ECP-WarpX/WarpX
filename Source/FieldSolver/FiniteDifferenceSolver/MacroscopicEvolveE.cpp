#include "FiniteDifferenceSolver.H"

#ifdef WARPX_DIM_RZ
    // currently works only for 3D
#else
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/FieldAccessorFunctors.H"
#endif
#include "MacroscopicProperties/MacroscopicProperties.H"
#include "Utils/CoarsenIO.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

using namespace amrex;

void FiniteDifferenceSolver::MacroscopicEvolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    amrex::Real const dt,
    std::unique_ptr<MacroscopicProperties> const& macroscopic_properties, int const lev)
{

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    amrex::ignore_unused(Efield, Bfield, Jfield, dt, macroscopic_properties, lev);
    amrex::Abort("currently macro E-push does not work for RZ");
#else
    if (m_do_nodal) {
        amrex::Abort(" macro E-push does not work for nodal ");

    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {

        if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff) {

            MacroscopicEvolveECartesian <CartesianYeeAlgorithm, LaxWendroffAlgo>
                       ( Efield, Bfield, Jfield, dt, macroscopic_properties, lev );

        }
        if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler) {

            MacroscopicEvolveECartesian <CartesianYeeAlgorithm, BackwardEulerAlgo>
                       ( Efield, Bfield, Jfield, dt, macroscopic_properties, lev );

        }

    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {

        // Note : EvolveE is the same for CKC and Yee.
        // In the templated Yee and CKC calls, the core operations for EvolveE is the same.
        if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff) {

            MacroscopicEvolveECartesian <CartesianCKCAlgorithm, LaxWendroffAlgo>
                       ( Efield, Bfield, Jfield, dt, macroscopic_properties, lev );

        } else if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler) {

            MacroscopicEvolveECartesian <CartesianCKCAlgorithm, BackwardEulerAlgo>
                       ( Efield, Bfield, Jfield, dt, macroscopic_properties, lev );

        }

    } else {
        amrex::Abort("MacroscopicEvolveE: Unknown algorithm");
    }
#endif

}


#ifndef WARPX_DIM_RZ

template<typename T_Algo, typename T_MacroAlgo>
void FiniteDifferenceSolver::MacroscopicEvolveECartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    amrex::Real const dt,
    std::unique_ptr<MacroscopicProperties> const& macroscopic_properties, int const lev)
{

    amrex::GpuArray<int, 3> const& Ex_stag = macroscopic_properties->Ex_IndexType;
    amrex::GpuArray<int, 3> const& Ey_stag = macroscopic_properties->Ey_IndexType;
    amrex::GpuArray<int, 3> const& Ez_stag = macroscopic_properties->Ez_IndexType;
    amrex::GpuArray<int, 3> const& Bx_stag = macroscopic_properties->Bx_IndexType;
    amrex::GpuArray<int, 3> const& By_stag = macroscopic_properties->By_IndexType;
    amrex::GpuArray<int, 3> const& Bz_stag = macroscopic_properties->Bz_IndexType;

    const auto getSigma = GetSigmaMacroparameter();
    const auto getEpsilon = GetEpsilonMacroparameter();
    const auto getMu = GetMuMacroparameter();
    auto &warpx = WarpX::GetInstance();
    const auto problo = warpx.Geom(lev).ProbLoArray();
    const auto dx = warpx.Geom(lev).CellSizeArray();

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
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
        Array4<Real> const& jx = Jfield[0]->array(mfi);
        Array4<Real> const& jy = Jfield[1]->array(mfi);
        Array4<Real> const& jz = Jfield[2]->array(mfi);

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        int const n_coefs_x = m_stencil_coefs_x.size();
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        int const n_coefs_y = m_stencil_coefs_y.size();
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        int const n_coefs_z = m_stencil_coefs_z.size();

        FieldAccessorMacroscopic<GetMuMacroparameter> const Hx(Bx, getMu, Bx_stag, problo, dx);
        FieldAccessorMacroscopic<GetMuMacroparameter> const Hy(By, getMu, By_stag, problo, dx);
        FieldAccessorMacroscopic<GetMuMacroparameter> const Hz(Bz, getMu, Bz_stag, problo, dx);

        // Extract tileboxes for which to loop
        Box const& tex  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        Box const& tey  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        Box const& tez  = mfi.tilebox(Efield[2]->ixType().toIntVect());
        // Loop over the cells and update the fields
        amrex::ParallelFor(tex, tey, tez,
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                amrex::Real x, y, z;
                WarpXUtilAlgo::getCellCoordinates (i, j, k, Ex_stag, problo, dx,
                                                   x, y, z );
                amrex::Real const sigma = getSigma(x, y, z);
                amrex::Real const epsilon = getEpsilon(x, y, z);
                amrex::Real alpha = T_MacroAlgo::alpha( sigma, epsilon, dt);
                amrex::Real beta = T_MacroAlgo::beta( sigma, epsilon, dt);
                Ex(i, j, k) = alpha * Ex(i, j, k)
                            + beta * ( - T_Algo::DownwardDz(Hy, coefs_z, n_coefs_z, i, j, k,0)
                                       + T_Algo::DownwardDy(Hz, coefs_y, n_coefs_y, i, j, k,0)
                                     ) - beta * jx(i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                amrex::Real x, y, z;
                WarpXUtilAlgo::getCellCoordinates (i, j, k, Ey_stag, problo, dx,
                                                   x, y, z );
                amrex::Real const sigma = getSigma(x, y, z);
                amrex::Real const epsilon = getEpsilon(x, y, z);
                amrex::Real alpha = T_MacroAlgo::alpha( sigma, epsilon, dt);
                amrex::Real beta = T_MacroAlgo::beta( sigma, epsilon, dt);

                Ey(i, j, k) = alpha * Ey(i, j, k)
                            + beta * ( - T_Algo::DownwardDx(Hz, coefs_x, n_coefs_x, i, j, k,0)
                                       + T_Algo::DownwardDz(Hx, coefs_z, n_coefs_z, i, j, k,0)
                                     ) - beta * jy(i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                amrex::Real x, y, z;
                WarpXUtilAlgo::getCellCoordinates (i, j, k, Ez_stag, problo, dx,
                                                   x, y, z );
                amrex::Real const sigma = getSigma(x, y, z);
                amrex::Real const epsilon = getEpsilon(x, y, z);
                amrex::Real alpha = T_MacroAlgo::alpha( sigma, epsilon, dt);
                amrex::Real beta = T_MacroAlgo::beta( sigma, epsilon, dt);

                Ez(i, j, k) = alpha * Ez(i, j, k)
                            + beta * ( - T_Algo::DownwardDy(Hx, coefs_y, n_coefs_y, i, j, k,0)
                                       + T_Algo::DownwardDx(Hy, coefs_x, n_coefs_x, i, j, k,0)
                                     ) - beta * jz(i, j, k);
            }
        );
    }
}

#endif // corresponds to ifndef WARPX_DIM_RZ
