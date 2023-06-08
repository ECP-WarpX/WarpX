/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FiniteDifferenceSolver.H"

#ifdef WARPX_DIM_RZ
#   include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#endif
#include "HybridPICModel/HybridPICModel.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>

using namespace amrex;

void FiniteDifferenceSolver::CalculateCurrentAmpere (
    std::array< std::unique_ptr<amrex::MultiFab>, 3>& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev )
{
   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
    if (m_fdtd_algo == ElectromagneticSolverAlgo::HybridPIC) {
#ifdef WARPX_DIM_RZ
        CalculateCurrentAmpereCylindrical <CylindricalYeeAlgorithm> (
            Jfield, Bfield, edge_lengths, lev
        );

#else
        CalculateCurrentAmpereCartesian <CartesianYeeAlgorithm> (
            Jfield, Bfield, edge_lengths, lev
        );

#endif
    } else {
        amrex::Abort(Utils::TextMsg::Err(
            "CalculateCurrentAmpere: Unknown algorithm choice."));
    }
}

// /**
//   * \brief Calculate electron current from Ampere's law without displacement
//   * current and the kinetically tracked ion currents i.e. J_e = curl B.
//   *
//   * \param[out] Jfield  vector of electron current MultiFabs at a given level
//   * \param[in] Bfield   vector of magnetic field MultiFabs at a given level
//   */
#ifdef WARPX_DIM_RZ
template<typename T_Algo>
void FiniteDifferenceSolver::CalculateCurrentAmpereCylindrical (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev
)
{
    amrex::ignore_unused(Jfield, Bfield, edge_lengths, lev);
    amrex::Abort(Utils::TextMsg::Err(
        "currently hybrid E-solve does not work for RZ"));
}
#else
template<typename T_Algo>
void FiniteDifferenceSolver::CalculateCurrentAmpereCartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev
)
{
    // for the profiler
    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif

    // reset Jfield
    Jfield[0]->setVal(0);
    Jfield[1]->setVal(0);
    Jfield[2]->setVal(0);

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Jfield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        // Extract field data for this grid/tile
        Array4<Real> const& Jx = Jfield[0]->array(mfi);
        Array4<Real> const& Jy = Jfield[1]->array(mfi);
        Array4<Real> const& Jz = Jfield[2]->array(mfi);
        Array4<Real const> const& Bx = Bfield[0]->const_array(mfi);
        Array4<Real const> const& By = Bfield[1]->const_array(mfi);
        Array4<Real const> const& Bz = Bfield[2]->const_array(mfi);

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

        Real const one_over_mu0 = 1._rt / PhysConst::mu0;

        // First calculate the total current using Ampere's law on the
        // same grid as the E-field
        amrex::ParallelFor(tjx, tjy, tjz,

            // Jx calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip if this cell is fully covered by embedded boundaries
                if (lx(i, j, k) <= 0) return;
#endif
                Jx(i, j, k) = one_over_mu0 * (
                    - T_Algo::DownwardDz(By, coefs_z, n_coefs_z, i, j, k)
                    + T_Algo::DownwardDy(Bz, coefs_y, n_coefs_y, i, j, k)
                );
            },

            // Jy calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip if this cell is fully covered by embedded boundaries
#ifdef WARPX_DIM_3D
                if (ly(i,j,k) <= 0) return;
#elif defined(WARPX_DIM_XZ)
                // In XZ Jy is associated with a mesh node, so we need to check if the mesh node is covered
                amrex::ignore_unused(ly);
                if (lx(i, j, k)<=0 || lx(i-1, j, k)<=0 || lz(i, j-1, k)<=0 || lz(i, j, k)<=0) return;
#endif
#endif
                Jy(i, j, k) = one_over_mu0 * (
                    - T_Algo::DownwardDx(Bz, coefs_x, n_coefs_x, i, j, k)
                    + T_Algo::DownwardDz(Bx, coefs_z, n_coefs_z, i, j, k)
                );
            },

            // Jz calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){
#ifdef AMREX_USE_EB
                // Skip if this cell is fully covered by embedded boundaries
                if (lz(i,j,k) <= 0) return;
#endif
                Jz(i, j, k) = one_over_mu0 * (
                    - T_Algo::DownwardDy(Bx, coefs_y, n_coefs_y, i, j, k)
                    + T_Algo::DownwardDx(By, coefs_x, n_coefs_x, i, j, k)
                );
            }
        );

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }
}
#endif


void FiniteDifferenceSolver::HybridPICSolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::unique_ptr<amrex::MultiFab> const& Pefield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, HybridPICModel const* hybrid_model,
    DtType a_dt_type )
{

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::grid_type == GridType::Staggered,
        "Ohm's law E-solve only works with a staggered (Yee) grid.");


   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
    if (m_fdtd_algo == ElectromagneticSolverAlgo::HybridPIC) {
#ifdef WARPX_DIM_RZ

        HybridPICSolveECylindrical <CylindricalYeeAlgorithm> (
            Efield, Jfield, Jifield, Bfield, rhofield, Pefield,
            edge_lengths, lev, hybrid_model, a_dt_type
        );

#else

        HybridPICSolveECartesian <CartesianYeeAlgorithm> (
            Efield, Jfield, Jifield, Bfield, rhofield, Pefield,
            edge_lengths, lev, hybrid_model, a_dt_type
        );

#endif
    } else {
        amrex::Abort(Utils::TextMsg::Err(
            "HybridSolveE: The hybrid-PIC electromagnetic solver algorithm must be used"));
    }
}

#ifdef WARPX_DIM_RZ
template<typename T_Algo>
void FiniteDifferenceSolver::HybridPICSolveECylindrical (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::unique_ptr<amrex::MultiFab> const& Pefield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, HybridPICModel const* hybrid_model,
    DtType a_dt_type )
{
#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif
    amrex::ignore_unused(
        Efield, Jfield, Jifield, Bfield, rhofield, Pefield, edge_lengths,
        lev, hybrid_model, a_dt_type
    );
    amrex::Abort(Utils::TextMsg::Err(
        "currently hybrid E-solve does not work for RZ"));
}

#else

template<typename T_Algo>
void FiniteDifferenceSolver::HybridPICSolveECartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jifield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::unique_ptr<amrex::MultiFab> const& Pefield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    int lev, HybridPICModel const* hybrid_model,
    DtType a_dt_type )
{
#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif

    // for the profiler
    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    using namespace ablastr::coarsen::sample;

    // get hybrid model parameters
    const auto eta = hybrid_model->m_eta;
    const auto rho_floor = hybrid_model->m_n_floor * PhysConst::q_e;

    // Index type required for interpolating fields from their respective
    // staggering to the Ex, Ey, Ez locations
    amrex::GpuArray<int, 3> const& Ex_stag = hybrid_model->Ex_IndexType;
    amrex::GpuArray<int, 3> const& Ey_stag = hybrid_model->Ey_IndexType;
    amrex::GpuArray<int, 3> const& Ez_stag = hybrid_model->Ez_IndexType;
    amrex::GpuArray<int, 3> const& Jx_stag = hybrid_model->Jx_IndexType;
    amrex::GpuArray<int, 3> const& Jy_stag = hybrid_model->Jy_IndexType;
    amrex::GpuArray<int, 3> const& Jz_stag = hybrid_model->Jz_IndexType;
    amrex::GpuArray<int, 3> const& Bx_stag = hybrid_model->Bx_IndexType;
    amrex::GpuArray<int, 3> const& By_stag = hybrid_model->By_IndexType;
    amrex::GpuArray<int, 3> const& Bz_stag = hybrid_model->Bz_IndexType;

    // Parameters for `interp` that maps from Yee to nodal mesh and back
    amrex::GpuArray<int, 3> const& nodal = {1, 1, 1};
    // The "coarsening is just 1 i.e. no coarsening"
    amrex::GpuArray<int, 3> const& coarsen = {1, 1, 1};

    // The E-field calculation is done in 2 steps:
    // 1) The J x B term is calculated on a nodal mesh in order to ensure
    //    energy conservation.
    // 2) The nodal E-field values are averaged onto the Yee grid and the
    //    electron pressure & resistivity terms are added (these terms are
    //    naturally located on the Yee grid).

    // Create a temporary multifab to hold the nodal E-field values
    // Note the multifab has 3 values for Ex, Ey and Ez which we can do here
    // since all three components will be calculated on the same grid.
    // Also note that enE_nodal_mf does not need to have any guard cells since
    // these values will be interpolated to the Yee mesh which is contained
    // by the nodal mesh.
    auto const& ba = convert(rhofield->boxArray(), IntVect::TheNodeVector());
    MultiFab enE_nodal_mf(ba, rhofield->DistributionMap(), 3, IntVect::TheZeroVector());

    // Loop through the grids, and over the tiles within each grid for the
    // initial, nodal calculation of E
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(enE_nodal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        Array4<Real> const& enE_nodal = enE_nodal_mf.array(mfi);
        Array4<Real const> const& Jx = Jfield[0]->const_array(mfi);
        Array4<Real const> const& Jy = Jfield[1]->const_array(mfi);
        Array4<Real const> const& Jz = Jfield[2]->const_array(mfi);
        Array4<Real const> const& Jix = Jifield[0]->const_array(mfi);
        Array4<Real const> const& Jiy = Jifield[1]->const_array(mfi);
        Array4<Real const> const& Jiz = Jifield[2]->const_array(mfi);
        Array4<Real const> const& Bx = Bfield[0]->const_array(mfi);
        Array4<Real const> const& By = Bfield[1]->const_array(mfi);
        Array4<Real const> const& Bz = Bfield[2]->const_array(mfi);

        // Loop over the cells and update the nodal E field
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // interpolate the total current to a nodal grid
            auto const jx_interp = Interp(Jx, Jx_stag, nodal, coarsen, i, j, k, 0);
            auto const jy_interp = Interp(Jy, Jy_stag, nodal, coarsen, i, j, k, 0);
            auto const jz_interp = Interp(Jz, Jz_stag, nodal, coarsen, i, j, k, 0);

            // interpolate the ion current to a nodal grid
            auto const jix_interp = Interp(Jix, Jx_stag, nodal, coarsen, i, j, k, 0);
            auto const jiy_interp = Interp(Jiy, Jy_stag, nodal, coarsen, i, j, k, 0);
            auto const jiz_interp = Interp(Jiz, Jz_stag, nodal, coarsen, i, j, k, 0);

            // interpolate the B field to a nodal grid
            auto const Bx_interp = Interp(Bx, Bx_stag, nodal, coarsen, i, j, k, 0);
            auto const By_interp = Interp(By, By_stag, nodal, coarsen, i, j, k, 0);
            auto const Bz_interp = Interp(Bz, Bz_stag, nodal, coarsen, i, j, k, 0);

            // calculate enE = (J - Ji) x B
            enE_nodal(i, j, k, 0) = (
                (jy_interp - jiy_interp) * Bz_interp
                - (jz_interp - jiz_interp) * By_interp
            );
            enE_nodal(i, j, k, 1) = (
                (jz_interp - jiz_interp) * Bx_interp
                - (jx_interp - jix_interp) * Bz_interp
            );
            enE_nodal(i, j, k, 2) = (
                (jx_interp - jix_interp) * By_interp
                - (jy_interp - jiy_interp) * Bx_interp
            );
        });

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }

    // Loop through the grids, and over the tiles within each grid again
    // for the Yee grid calculation of the E field
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
        Array4<Real const> const& Jx = Jfield[0]->const_array(mfi);
        Array4<Real const> const& Jy = Jfield[1]->const_array(mfi);
        Array4<Real const> const& Jz = Jfield[2]->const_array(mfi);
        Array4<Real const> const& enE = enE_nodal_mf.const_array(mfi);
        Array4<Real const> const& rho = rhofield->const_array(mfi);
        Array4<Real> const& Pe = Pefield->array(mfi);

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
                // Interpolate to get the appropriate charge density in space
                Real rho_val = Interp(rho, nodal, Ex_stag, coarsen, i, j, k, 0);

                // safety condition since we divide by rho_val later
                if (rho_val < rho_floor) rho_val = rho_floor;

                // Get the gradient of the electron pressure
                auto grad_Pe = T_Algo::UpwardDx(Pe, coefs_x, n_coefs_x, i, j, k);

                // interpolate the nodal neE values to the Yee grid
                auto enE_x = Interp(enE, nodal, Ex_stag, coarsen, i, j, k, 0);

                Ex(i, j, k) = (enE_x - grad_Pe) / rho_val;

                // Add resistivity only if E field value is used to update B
                if (a_dt_type != DtType::Full) Ex(i, j, k) += eta(rho_val) * Jx(i, j, k);
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
                // Interpolate to get the appropriate charge density in space
                Real rho_val = Interp(rho, nodal, Ey_stag, coarsen, i, j, k, 0);

                // safety condition since we divide by rho_val later
                if (rho_val < rho_floor) rho_val = rho_floor;

                // Get the gradient of the electron pressure
                auto grad_Pe = T_Algo::UpwardDy(Pe, coefs_y, n_coefs_y, i, j, k);

                // interpolate the nodal neE values to the Yee grid
                auto enE_y = Interp(enE, nodal, Ey_stag, coarsen, i, j, k, 1);

                Ey(i, j, k) = (enE_y - grad_Pe) / rho_val;

                // Add resistivity only if E field value is used to update B
                if (a_dt_type != DtType::Full) Ey(i, j, k) += eta(rho_val) * Jy(i, j, k);
            },

            // Ez calculation
            [=] AMREX_GPU_DEVICE (int i, int j, int k){

#ifdef AMREX_USE_EB
                // Skip field solve if this cell is fully covered by embedded boundaries
                if (lz(i,j,k) <= 0) return;
#endif
                // Interpolate to get the appropriate charge density in space
                Real rho_val = Interp(rho, nodal, Ez_stag, coarsen, i, j, k, 0);

                // safety condition since we divide by rho_val later
                if (rho_val < rho_floor) rho_val = rho_floor;

                // Get the gradient of the electron pressure
                auto grad_Pe = T_Algo::UpwardDz(Pe, coefs_z, n_coefs_z, i, j, k);

                // interpolate the nodal neE values to the Yee grid
                auto enE_z = Interp(enE, nodal, Ez_stag, coarsen, i, j, k, 2);

                Ez(i, j, k) = (enE_z - grad_Pe) / rho_val;

                // Add resistivity only if E field value is used to update B
                if (a_dt_type != DtType::Full) Ez(i, j, k) += eta(rho_val) * Jz(i, j, k);
            }
        );

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }
}
#endif
