#include "FiniteDifferenceSolver.H"

#ifdef WARPX_DIM_RZ
    // currently works only for 3D
#else
#   include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#   include "FiniteDifferenceAlgorithms/FieldAccessorFunctors.H"
#endif
#include "MacroscopicProperties/MacroscopicProperties.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>

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

void FiniteDifferenceSolver::MacroscopicEvolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Hfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Mfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    amrex::Real const dt,
    std::unique_ptr<MacroscopicProperties> const& macroscopic_properties)
{

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    amrex::ignore_unused(Efield, Bfield, Jfield, edge_lengths, dt, macroscopic_properties);

    WARPX_ABORT_WITH_MESSAGE("currently macro E-push does not work for RZ");
#else
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_grid_type != GridType::Collocated, "Macroscopic E field solver does not work on collocated grids");


    if (m_fdtd_algo == ElectromagneticSolverAlgo::Yee) {

        if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff) {

            MacroscopicEvolveECartesian <CartesianYeeAlgorithm, LaxWendroffAlgo>
                       ( Efield, Bfield, Hfield, Mfield, Jfield, edge_lengths, dt, macroscopic_properties);

        }
        if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler) {

            MacroscopicEvolveECartesian <CartesianYeeAlgorithm, BackwardEulerAlgo>
                       ( Efield, Bfield, Hfield, Mfield, Jfield, edge_lengths, dt, macroscopic_properties);

        }

    } else if (m_fdtd_algo == ElectromagneticSolverAlgo::CKC) {

        // Note : EvolveE is the same for CKC and Yee.
        // In the templated Yee and CKC calls, the core operations for EvolveE is the same.
        if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff) {

            MacroscopicEvolveECartesian <CartesianCKCAlgorithm, LaxWendroffAlgo>
                       ( Efield, Bfield, Hfield, Mfield, Jfield, edge_lengths, dt, macroscopic_properties);

        } else if (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler) {

            MacroscopicEvolveECartesian <CartesianCKCAlgorithm, BackwardEulerAlgo>
                       ( Efield, Bfield, Hfield, Mfield, Jfield, edge_lengths, dt, macroscopic_properties);

        }

    } else {
        WARPX_ABORT_WITH_MESSAGE(
            "MacroscopicEvolveE: Unknown algorithm");
    }
#endif

}


#ifndef WARPX_DIM_RZ

template<typename T_Algo, typename T_MacroAlgo>
void FiniteDifferenceSolver::MacroscopicEvolveECartesian (
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Hfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Mfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
    amrex::Real const dt,
    std::unique_ptr<MacroscopicProperties> const& macroscopic_properties)
{
    using Real = amrex::Real;

    // using T_Algo = CartesianYeeAlgorithm; //! debug variable
    // using T_MacroAlgo = BackwardEulerAlgo; //! debug variable

#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif

    amrex::MultiFab const& sigma_mf = macroscopic_properties->getsigma_mf();
    amrex::MultiFab const& epsilon_mf = macroscopic_properties->getepsilon_mf();
    amrex::MultiFab const& mu_mf = macroscopic_properties->getmu_mf();
    amrex::iMultiFab const& mag_mat_id_mf = macroscopic_properties->get_magnetic_material_id_mf();

    // Index type required for calling ablastr::coarsen::sample::Interp to interpolate macroscopic
    // properties from their respective staggering to the Ex, Ey, Ez locations
    amrex::GpuArray<int, 3> const& sigma_stag = macroscopic_properties->sigma_IndexType;
    amrex::GpuArray<int, 3> const& epsilon_stag = macroscopic_properties->epsilon_IndexType;
    amrex::GpuArray<int, 3> const& macro_cr     = macroscopic_properties->macro_cr_ratio;
    amrex::GpuArray<int, 3> const& Ex_stag = macroscopic_properties->Ex_IndexType;
    amrex::GpuArray<int, 3> const& Ey_stag = macroscopic_properties->Ey_IndexType;
    amrex::GpuArray<int, 3> const& Ez_stag = macroscopic_properties->Ez_IndexType;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(*Efield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        amrex::Array4<      Real> const& Ex = Efield[0]->array(mfi);
        amrex::Array4<      Real> const& Ey = Efield[1]->array(mfi);
        amrex::Array4<      Real> const& Ez = Efield[2]->array(mfi);
        amrex::Array4<const Real> const& Bx = Bfield[0]->const_array(mfi);
        amrex::Array4<const Real> const& By = Bfield[1]->const_array(mfi);
        amrex::Array4<const Real> const& Bz = Bfield[2]->const_array(mfi);
        amrex::Array4<const Real> const& jx = Jfield[0]->const_array(mfi);
        amrex::Array4<const Real> const& jy = Jfield[1]->const_array(mfi);
        amrex::Array4<const Real> const& jz = Jfield[2]->const_array(mfi);

#ifdef AMREX_USE_EB
        amrex::Array4<const Real> const& lx = edge_lengths[0]->const_array(mfi);
        amrex::Array4<const Real> const& ly = edge_lengths[1]->const_array(mfi);
        amrex::Array4<const Real> const& lz = edge_lengths[2]->const_array(mfi);
#endif

        // material prop //
        amrex::Array4<const Real> const& sigma_arr = sigma_mf.const_array(mfi);
        amrex::Array4<const Real> const& eps_arr = epsilon_mf.const_array(mfi);
        amrex::Array4<const Real> const& mu_arr = mu_mf.const_array(mfi);

        // Extract stencil coefficients
        Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
        auto const n_coefs_x = static_cast<int>(m_stencil_coefs_x.size());
        Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
        auto const n_coefs_y = static_cast<int>(m_stencil_coefs_y.size());
        Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
        auto const n_coefs_z = static_cast<int>(m_stencil_coefs_z.size());

        auto Hx = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) { return Real(0.0); };
        auto Hy = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) { return Real(0.0); };
        auto Hz = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) { return Real(0.0); };

        if (macroscopic_properties->is_magnetic_material_present()) {
            amrex::Array4<Real> const& H_x_arr = Hfield[0]->array(mfi);
            amrex::Array4<Real> const& H_y_arr = Hfield[1]->array(mfi);
            amrex::Array4<Real> const& H_z_arr = Hfield[2]->array(mfi);
            amrex::Array4<Real> const& M_x_arr = Mfield[0]->array(mfi);
            amrex::Array4<Real> const& M_y_arr = Mfield[1]->array(mfi);
            amrex::Array4<Real> const& M_z_arr = Mfield[2]->array(mfi);

            amrex::Array4<const int> const& mag_mat_id = mag_mat_id_mf.const_array(mfi);
            amrex::Vector<const MagneticMaterial> const& mag_mat_vector = macroscopic_properties->m_magnetic_materials;

            amrex::Box const& box_cc = mfi.tilebox(amrex::IntVect::TheCellVector());
            amrex::ParallelFor(box_cc, 1,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    if (mag_mat_id(i,j,k) >= 0) {
                        const Real B_next_x = Real(0.5)*(Bx(i,j,k,n) + Bx(i+1,j,k,n));
                        const Real B_next_y = Real(0.5)*(By(i,j,k,n) + By(i,j+1,k,n));
                        const Real B_next_z = Real(0.5)*(Bz(i,j,k,n) + Bz(i,j,k+1,n));
                        mag_mat_vector[mag_mat_id(i,j,k)].UpdateHandM(H_x_arr(i,j,k,n), H_y_arr(i,j,k,n), H_z_arr(i,j,k,n),
                                                                      M_x_arr(i,j,k,n), M_y_arr(i,j,k,n), M_z_arr(i,j,k,n),
                                                                      B_next_x        , B_next_y        , B_next_z        );
                    }
                }
            );
            Hx = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                return Real(0.5*((mag_mat_id(i  ,j,k) >= 0 ?  H_x_arr(i  ,j,k,n) : Bx(i,j,k,n)/mu_arr(i  ,j,k))
                                +(mag_mat_id(i-1,j,k) >= 0 ?  H_x_arr(i-1,j,k,n) : Bx(i,j,k,n)/mu_arr(i-1,j,k))));
            };
            Hy = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                return Real(0.5*((mag_mat_id(i,j  ,k) >= 0 ?  H_y_arr(i,j  ,k,n) : By(i,j,k,n)/mu_arr(i,j  ,k))
                                +(mag_mat_id(i,j-1,k) >= 0 ?  H_y_arr(i,j-1,k,n) : By(i,j,k,n)/mu_arr(i,j-1,k))));
            };
            Hz = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                return Real(0.5*((mag_mat_id(i,j,k  ) >= 0 ?  H_z_arr(i,j,k  ,n) : Bz(i,j,k,n)/mu_arr(i,j,k  ))
                                +(mag_mat_id(i,j,k-1) >= 0 ?  H_z_arr(i,j,k-1,n) : Bz(i,j,k,n)/mu_arr(i,j,k-1))));
            };
        } else {
            Hx = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                return (Real(0.5)/mu_arr(i,j,k) + Real(0.5)/mu_arr(i-1,j,k))*Bx(i,j,k,n);
            };
            Hy = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                return (Real(0.5)/mu_arr(i,j,k) + Real(0.5)/mu_arr(i,j-1,k))*By(i,j,k,n);
            };
            Hz = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                return (Real(0.5)/mu_arr(i,j,k) + Real(0.5)/mu_arr(i,j,k-1))*Bz(i,j,k,n);
            };
        }

        // Extract tileboxes for which to loop
        amrex::Box const& tex  = mfi.tilebox(Efield[0]->ixType().toIntVect());
        amrex::Box const& tey  = mfi.tilebox(Efield[1]->ixType().toIntVect());
        amrex::Box const& tez  = mfi.tilebox(Efield[2]->ixType().toIntVect());
        // starting component to interpolate macro properties to Ex, Ey, Ez locations
        const int scomp = 0;
        // Loop over the cells and update the fields
        amrex::ParallelFor(tex, tey, tez,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef AMREX_USE_EB
                // Skip field push if this cell is fully covered by embedded boundaries
                if (lx(i, j, k) <= 0) return;
#endif
                // Interpolate conductivity, sigma, to Ex position on the grid
                amrex::Real const sigma_interp = ablastr::coarsen::sample::Interp(sigma_arr, sigma_stag,
                                                                                  Ex_stag, macro_cr, i, j, k, scomp);
                // Interpolated permittivity, epsilon, to Ex position on the grid
                amrex::Real const epsilon_interp = ablastr::coarsen::sample::Interp(eps_arr, epsilon_stag,
                                                                                    Ex_stag, macro_cr, i, j, k, scomp);
                const amrex::Real alpha = T_MacroAlgo::alpha( sigma_interp, epsilon_interp, dt);
                const amrex::Real beta = T_MacroAlgo::beta( sigma_interp, epsilon_interp, dt);
                Ex(i, j, k) = alpha * Ex(i, j, k)
                            + beta * ( - T_Algo::DownwardDz(Hy, coefs_z, n_coefs_z, i, j, k, 0)
                                       + T_Algo::DownwardDy(Hz, coefs_y, n_coefs_y, i, j, k, 0)
                                       - jx(i, j, k) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef AMREX_USE_EB
#ifdef WARPX_DIM_3D
                if (ly(i,j,k) <= 0) return;
#elif defined(WARPX_DIM_XZ)
                //In XZ Ey is associated with a mesh node, so we need to check if the mesh node is covered
                amrex::ignore_unused(ly);
                if (lx(i, j, k)<=0 || lx(i-1, j, k)<=0 || lz(i, j, k)<=0 || lz(i, j-1, k)<=0) return;
#endif
#endif
                // Interpolate conductivity, sigma, to Ey position on the grid
                amrex::Real const sigma_interp = ablastr::coarsen::sample::Interp(sigma_arr, sigma_stag,
                                                                                  Ey_stag, macro_cr, i, j, k, scomp);
                // Interpolated permittivity, epsilon, to Ey position on the grid
                amrex::Real const epsilon_interp = ablastr::coarsen::sample::Interp(eps_arr, epsilon_stag,
                                                                                    Ey_stag, macro_cr, i, j, k, scomp);
                const amrex::Real alpha = T_MacroAlgo::alpha( sigma_interp, epsilon_interp, dt);
                const amrex::Real beta = T_MacroAlgo::beta( sigma_interp, epsilon_interp, dt);

                Ey(i, j, k) = alpha * Ey(i, j, k)
                            + beta * ( - T_Algo::DownwardDx(Hz, coefs_x, n_coefs_x, i, j, k, 0)
                                       + T_Algo::DownwardDz(Hx, coefs_z, n_coefs_z, i, j, k, 0)
                                       - jy(i, j, k) );
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef AMREX_USE_EB
                // Skip field push if this cell is fully covered by embedded boundaries
                if (lz(i,j,k) <= 0) return;
#endif
                // Interpolate conductivity, sigma, to Ez position on the grid
                amrex::Real const sigma_interp = ablastr::coarsen::sample::Interp(sigma_arr, sigma_stag,
                                                                                  Ez_stag, macro_cr, i, j, k, scomp);
                // Interpolated permittivity, epsilon, to Ez position on the grid
                amrex::Real const epsilon_interp = ablastr::coarsen::sample::Interp(eps_arr, epsilon_stag,
                                                                                    Ez_stag, macro_cr, i, j, k, scomp);
                const amrex::Real alpha = T_MacroAlgo::alpha( sigma_interp, epsilon_interp, dt);
                const amrex::Real beta = T_MacroAlgo::beta( sigma_interp, epsilon_interp, dt);

                Ez(i, j, k) = alpha * Ez(i, j, k)
                            + beta * ( - T_Algo::DownwardDy(Hx, coefs_y, n_coefs_y, i, j, k, 0)
                                       + T_Algo::DownwardDx(Hy, coefs_x, n_coefs_x, i, j, k, 0)
                                       - jz(i, j, k) );
            }
        );
    }
}

#endif // corresponds to ifndef WARPX_DIM_RZ
