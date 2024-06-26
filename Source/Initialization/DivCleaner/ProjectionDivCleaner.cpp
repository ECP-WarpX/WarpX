 /* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: S. Eric Clark (Helion Energy, Inc.)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ProjectionDivCleaner.H"

#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_LO_BCTYPES.H>

#include <WarpX.H>
#if defined WARPX_DIM_RZ
#include <FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H>
#else
#include <FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H>
#endif

#include <ablastr/utils/Communication.H>

using namespace amrex;

namespace warpx::initialization {

ProjectionDivCleaner::ProjectionDivCleaner()
{
    // Initialize tolerance based on field precision
    if constexpr (std::is_same<Real, float>::value) {
        // Error out of divergence cleaner
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
            "Single Precision Divergence Cleaner has convergence problems. "
            "Please compile with WarpX_PRECISION=DOUBLE."
        );

        m_rtol = 1e-6;
        m_atol = 0.0;
    }
    else {
        m_rtol = 1e-12;
        m_atol = 0.0;
    }

    auto& warpx = WarpX::GetInstance();
    m_levels = warpx.finestLevel() + 1;

    m_solution.resize(m_levels);
    m_source.resize(m_levels);

    const int ncomps = WarpX::ncomps;
    auto const& ng = warpx.getFieldPointer(warpx::fields::FieldType::Bfield_aux, 0, 0)->nGrowVect();


    std::array<amrex::Real,3> const& cell_size = WarpX::CellSize(0);

    for (int lev = 0; lev < m_levels; ++lev)
    {
        // Default BoxArray and DistributionMap for initializing the output MultiFab, m_mf_output.
        const amrex::BoxArray& ba = warpx.boxArray(lev);
        const amrex::DistributionMapping& dmap = warpx.DistributionMap(lev);

        m_solution[lev].reset();
        m_source[lev].reset();

        const auto tag1 = amrex::MFInfo().SetTag("div_cleaner_solution");
        m_solution[lev] = std::make_unique<MultiFab>(amrex::convert(ba, IntVect::TheCellVector()),
            dmap, ncomps, ng, tag1);
        const auto tag2 = amrex::MFInfo().SetTag("div_cleaner_source");
        m_source[lev] = std::make_unique<MultiFab>(amrex::convert(ba, IntVect::TheCellVector()),
            dmap, ncomps, ng, tag2);

        m_solution[lev]->setVal(0.0);
        m_source[lev]->setVal(0.0);
    }

    m_h_stencil_coefs_x.resize(1);
    m_h_stencil_coefs_x[0] = 1._rt/cell_size[0];
    m_h_stencil_coefs_y.resize(1);
    m_h_stencil_coefs_y[0] = 1._rt/cell_size[1];
    m_h_stencil_coefs_z.resize(1);
    m_h_stencil_coefs_z[0] = 1._rt/cell_size[2];

    m_stencil_coefs_x.resize(m_h_stencil_coefs_x.size());
    m_stencil_coefs_y.resize(m_h_stencil_coefs_y.size());
    m_stencil_coefs_z.resize(m_h_stencil_coefs_z.size());

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                          m_h_stencil_coefs_x.begin(), m_h_stencil_coefs_x.end(),
                          m_stencil_coefs_x.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                          m_h_stencil_coefs_y.begin(), m_h_stencil_coefs_y.end(),
                          m_stencil_coefs_y.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                          m_h_stencil_coefs_z.begin(), m_h_stencil_coefs_z.end(),
                          m_stencil_coefs_z.begin());
    amrex::Gpu::synchronize();
}

void
ProjectionDivCleaner::solve ()
{
    // Get WarpX object
    auto & warpx = WarpX::GetInstance();

    const auto& ba = warpx.boxArray();
    const auto& dmap = warpx.DistributionMap();
    const auto& geom = warpx.Geom();

    // Pull boundary conditions from WarpX class
    // bogus values are overwritten.
    amrex::Array<LinOpBCType,AMREX_SPACEDIM> lobc({AMREX_D_DECL(LinOpBCType::bogus,
                                                                LinOpBCType::bogus,
                                                                LinOpBCType::bogus)});
    amrex::Array<LinOpBCType,AMREX_SPACEDIM> hibc({AMREX_D_DECL(LinOpBCType::bogus,
                                                                LinOpBCType::bogus,
                                                                LinOpBCType::bogus)});

#ifdef WARPX_DIM_RZ
    if (geom[0].ProbLo(0) == 0){
        lobc[0] = LinOpBCType::Neumann;

        // handle the r_max boundary explicitly
        if (WarpX::field_boundary_hi[0] == FieldBoundaryType::PEC) {
            hibc[0] = LinOpBCType::Dirichlet;
        }
        else if (WarpX::field_boundary_hi[0] == FieldBoundaryType::Neumann) {
            hibc[0] = LinOpBCType::Neumann;
        }
        else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                "Field boundary condition at the outer radius must be either PEC or neumann "
                "when using the Projection Divergence Cleaner"
            );
        }
    }
    const int dim_start = 1;
#else
    const int dim_start = 0;
#endif
    for (int idim=dim_start; idim<AMREX_SPACEDIM; idim++){
        if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::Periodic
             && WarpX::field_boundary_hi[idim] == FieldBoundaryType::Periodic ) {
            lobc[idim] = LinOpBCType::Periodic;
            hibc[idim] = LinOpBCType::Periodic;
        }
        else {
            if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC ) {
                lobc[idim] = LinOpBCType::Dirichlet;
            }
            else if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::Neumann ) {
                lobc[idim] = LinOpBCType::Neumann;
            }
            else {
                WARPX_ABORT_WITH_MESSAGE(
                    "Field boundary conditions have to be either periodic, PEC or neumann "
                    "when using the Projection based Divergence Cleaner,  but they are " + GetFieldBCTypeString(WarpX::field_boundary_lo[idim])
                );
            }

            if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC ) {
                hibc[idim] = LinOpBCType::Dirichlet;
            }
            else if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::Neumann ) {
                hibc[idim] = LinOpBCType::Neumann;
            }
            else {
                WARPX_ABORT_WITH_MESSAGE(
                    "Field boundary conditions have to be either periodic, PEC or neumann "
                    "when using the electrostatic Multigrid solver,  but they are " + GetFieldBCTypeString(WarpX::field_boundary_hi[idim])
                );
            }
        }

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            (WarpX::field_boundary_lo[idim] != FieldBoundaryType::Open &&
            WarpX::field_boundary_hi[idim] != FieldBoundaryType::Open &&
            WarpX::field_boundary_lo[idim] != FieldBoundaryType::PML &&
            WarpX::field_boundary_hi[idim] != FieldBoundaryType::PML) ,
            "Open and PML field boundary conditions only work with "
            "warpx.poisson_solver = fft."
        );
    }

    LPInfo info;
    info.setAgglomeration(m_agglomeration);
    info.setConsolidation(m_consolidation);
    info.setMaxCoarseningLevel(m_max_coarsening_level);

    for (int ilev = 0; ilev < m_levels; ++ilev)
    {
        MLPoisson mlpoisson({geom[ilev]}, {ba[ilev]}, {dmap[ilev]}, info);

        mlpoisson.setMaxOrder(m_linop_maxorder);
        mlpoisson.setDomainBC(lobc, hibc);

        if (ilev > 0) {
            mlpoisson.setCoarseFineBC(m_solution[ilev-1].get(), m_ref_ratio);
        }

        m_solution[ilev]->setVal(0.);

        mlpoisson.setLevelBC(ilev, m_solution[ilev].get());

        MLMG mlmg(mlpoisson);
        mlmg.setMaxIter(m_max_iter);
        mlmg.setMaxFmgIter(m_max_fmg_iter);
        mlmg.setVerbose(m_verbose);
        mlmg.setBottomVerbose(m_bottom_verbose);
        mlmg.solve({m_solution[ilev].get()}, {m_source[ilev].get()}, m_rtol, m_atol);

        // Synchronize the ghost cells, do halo exchange
        ablastr::utils::communication::FillBoundary(*m_solution[ilev],
                                                m_solution[ilev]->nGrowVect(),
                                                WarpX::do_single_precision_comms,
                                                geom[ilev].periodicity());
    }
}

void
ProjectionDivCleaner::setSourceFromBfield ()
{
    // Get WarpX object
    auto & warpx = WarpX::GetInstance();
    const auto& geom = warpx.Geom();

    // This function will compute -divB and store it in the source multifab
    for (int ilev = 0; ilev < m_levels; ++ilev)
    {
        // Zero out source multifab
        m_source[ilev]->setVal(0.0);

        WarpX::ComputeDivB(
            *m_source[ilev],
            0,
            warpx.getFieldPointerArray(warpx::fields::FieldType::Bfield_aux, ilev),
            WarpX::CellSize(0)
            );

        m_source[ilev]->mult(-1._rt);

        // Synchronize the ghost cells, do halo exchange
        ablastr::utils::communication::FillBoundary(*m_source[ilev],
                                                m_source[ilev]->nGrowVect(),
                                                WarpX::do_single_precision_comms,
                                                geom[ilev].periodicity());
    }
}

void
ProjectionDivCleaner::correctBfield ()
{
    // Get WarpX object
    auto & warpx = WarpX::GetInstance();
    const auto& geom = warpx.Geom();

    // This function computes the gradient of the solution and subtracts out divB component from B
    for (int ilev = 0; ilev < m_levels; ++ilev)
    {
        // Grab B-field multifabs at this level
        amrex::MultiFab* Bx = warpx.getFieldPointer(warpx::fields::FieldType::Bfield_aux, ilev, 0);
        amrex::MultiFab* By = warpx.getFieldPointer(warpx::fields::FieldType::Bfield_aux, ilev, 1);
        amrex::MultiFab* Bz = warpx.getFieldPointer(warpx::fields::FieldType::Bfield_aux, ilev, 2);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_solution[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Grab references to B field arrays for this grid/tile
            amrex::Array4<Real> const& Bx_arr = Bx->array(mfi);
#if !defined(WARPX_DIM_RZ)
            amrex::Array4<Real> const& By_arr = By->array(mfi);
#endif
            amrex::Array4<Real> const& Bz_arr = Bz->array(mfi);

            // Extract stencil coefficients
            Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
            auto const n_coefs_x = static_cast<int>(m_stencil_coefs_x.size());
#if !defined(WARPX_DIM_RZ)
            Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
            auto const n_coefs_y = static_cast<int>(m_stencil_coefs_y.size());
#endif
            Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
            auto const n_coefs_z = static_cast<int>(m_stencil_coefs_z.size());

            const Box& tbx = mfi.tilebox(Bx->ixType().toIntVect());
#if !defined(WARPX_DIM_RZ)
            const Box& tby = mfi.tilebox(By->ixType().toIntVect());
#endif
            const Box& tbz = mfi.tilebox(Bz->ixType().toIntVect());

            amrex::Array4<Real> const& sol_arr = m_solution[ilev]->array(mfi);
#if defined(WARPX_DIM_RZ)
            amrex::ParallelFor(tbx, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
            {
                Bx_arr(i,j,0) += CylindricalYeeAlgorithm::DownwardDr(sol_arr, coefs_x, n_coefs_x, i, j, 0, 0);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
            {
                Bz_arr(i,j,0) += CylindricalYeeAlgorithm::DownwardDz(sol_arr, coefs_z, n_coefs_z, i, j, 0, 0);
            });
#else
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Bx_arr(i,j,k) += CartesianYeeAlgorithm::DownwardDx(sol_arr, coefs_x, n_coefs_x, i, j, k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                By_arr(i,j,k) += CartesianYeeAlgorithm::DownwardDy(sol_arr, coefs_y, n_coefs_y, i, j, k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Bz_arr(i,j,k) += CartesianYeeAlgorithm::DownwardDz(sol_arr, coefs_z, n_coefs_z, i, j, k);
            });
#endif
        }
        // Synchronize the ghost cells, do halo exchange
        ablastr::utils::communication::FillBoundary(*Bx,
                                                    Bx->nGrowVect(),
                                                    WarpX::do_single_precision_comms,
                                                    geom[ilev].periodicity());
        ablastr::utils::communication::FillBoundary(*By,
                                                    By->nGrowVect(),
                                                    WarpX::do_single_precision_comms,
                                                    geom[ilev].periodicity());
        ablastr::utils::communication::FillBoundary(*Bz,
                                                    Bz->nGrowVect(),
                                                    WarpX::do_single_precision_comms,
                                                    geom[ilev].periodicity());
        amrex::Gpu::synchronize();
    }
}

} // namespace warpx::initialization

void
WarpX::ProjectionCleanDivB() {
#if defined(WARPX_DIM_RZ)
    ablastr::warn_manager::WMRecordWarning("Projection Div Cleaner",
        "WarpX is running in RZ mode, so divB not cleaned. Interpolation may lead to non-zero B field divergence.",
        ablastr::warn_manager::WarnPriority::low);
#else
    if constexpr (!std::is_same<Real, double>::value) {
        ablastr::warn_manager::WMRecordWarning("Projection Div Cleaner",
            "Field Precision is SINGLE, so divB not cleaned. Interpolation may lead to non-zero B field divergence.",
            ablastr::warn_manager::WarnPriority::low);
    } else if (grid_type == GridType::Collocated) {
        ablastr::warn_manager::WMRecordWarning("Projection Div Cleaner",
            "Grid Type is collocated, so divB not cleaned. Interpolation may lead to non-zero B field divergence.",
            ablastr::warn_manager::WarnPriority::low);
    } else if ( WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::Yee
            ||  WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::HybridPIC
            ||  ( (WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame
                || WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic)
                && WarpX::poisson_solver_id == PoissonSolverAlgo::Multigrid)) {
        // Build Object, run, then delete
        warpx::initialization::ProjectionDivCleaner dc;
        dc.setSourceFromBfield();
        dc.solve();
        dc.correctBfield();
    } else {
        ablastr::warn_manager::WMRecordWarning("Projection Div Cleaner",
            "Only Yee, HybridPIC, and MLMG based static Labframe solvers are currently supported, so divB not cleaned. Interpolation may lead to non-zero B field divergence.",
            ablastr::warn_manager::WarnPriority::low);
    }
#endif
}
