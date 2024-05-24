#include "ProjectionDivCleaner.H"

#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>

#include <WarpX.H>
#include <FieldSolver/FiniteDifferenceSolver/FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H>

#include <ablastr/utils/Communication.H>

using namespace amrex;

ProjectionDivCleaner::ProjectionDivCleaner()
{
    // Initialize tolerance based on field precision
    if constexpr (std::is_same<Real, float>::value) {
        // Error out of divergence cleaner
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
            "Single Precision Divergence Cleaner has convergence problems. "
            "Please compile with WarpX_PRECISION=DOUBLE"
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

    int ncomps = WarpX::ncomps;
    // auto Bx_point = warpx.getFieldPointer(warpx::fields::FieldType::Bfield_aux, 0, 0);
    auto const& ng = IntVect( AMREX_D_DECL(2,2,2) );


    std::array<amrex::Real,3> const& cell_size = warpx.CellSize(0);

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

    LPInfo info;
    info.setAgglomeration(m_agglomeration);
    info.setConsolidation(m_consolidation);
    info.setMaxCoarseningLevel(m_max_coarsening_level);

    for (int ilev = 0; ilev < m_levels; ++ilev)
    {
        MLPoisson mlpoisson({geom[ilev]}, {ba[ilev]}, {dmap[ilev]}, info);

        mlpoisson.setMaxOrder(m_linop_maxorder);

        // This is a 3d problem with Dirichlet BC
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                            LinOpBCType::Periodic,
                                            LinOpBCType::Neumann)},
                              {AMREX_D_DECL(LinOpBCType::Periodic,
                                            LinOpBCType::Periodic,
                                            LinOpBCType::Neumann)});

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
        ablastr::utils::communication::FillBoundary(*(m_solution[ilev].get()),
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

    // This function will compute divB and add to the source multifab
    for (int ilev = 0; ilev < m_levels; ++ilev)
    {
        // Zero out source multifab
        m_source[ilev]->setVal(0.0);

        // Grab B-field multifabs at this level
        std::array<const amrex::MultiFab* const, 3> Bfield =
            warpx.getFieldPointerArray(warpx::fields::FieldType::Bfield_aux, ilev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*(m_source[ilev].get()), TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Grab references to B field arrays for this grid/tile
            amrex::Array4<Real const> const& Bx_arr = Bfield[0]->const_array(mfi);
            amrex::Array4<Real const> const& By_arr = Bfield[1]->const_array(mfi);
            amrex::Array4<Real const> const& Bz_arr = Bfield[2]->const_array(mfi);

            // Extract stencil coefficients
            Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
            auto const n_coefs_x = static_cast<int>(m_stencil_coefs_x.size());
            Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
            auto const n_coefs_y = static_cast<int>(m_stencil_coefs_y.size());
            Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
            auto const n_coefs_z = static_cast<int>(m_stencil_coefs_z.size());

            amrex::Array4<Real> const& src_arr = m_source[ilev].get()->array(mfi);

            amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                src_arr(i,j,k) -= CartesianYeeAlgorithm::DownwardDx(Bx_arr, coefs_x, n_coefs_x, i+1, j, k);
                src_arr(i,j,k) -= CartesianYeeAlgorithm::DownwardDy(By_arr, coefs_y, n_coefs_y, i, j+1, k);
                src_arr(i,j,k) -= CartesianYeeAlgorithm::DownwardDz(Bz_arr, coefs_z, n_coefs_z, i, j, k+1);
            });
        }

        // Synchronize the ghost cells, do halo exchange
        ablastr::utils::communication::FillBoundary(*(m_source[ilev].get()),
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

    // This function will compute divB and add to the source multifab
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
            amrex::Array4<Real> const& By_arr = By->array(mfi);
            amrex::Array4<Real> const& Bz_arr = Bz->array(mfi);

            // Extract stencil coefficients
            Real const * const AMREX_RESTRICT coefs_x = m_stencil_coefs_x.dataPtr();
            auto const n_coefs_x = static_cast<int>(m_stencil_coefs_x.size());
            Real const * const AMREX_RESTRICT coefs_y = m_stencil_coefs_y.dataPtr();
            auto const n_coefs_y = static_cast<int>(m_stencil_coefs_y.size());
            Real const * const AMREX_RESTRICT coefs_z = m_stencil_coefs_z.dataPtr();
            auto const n_coefs_z = static_cast<int>(m_stencil_coefs_z.size());

            const Box& tbx = mfi.tilebox(Bx->ixType().toIntVect());
            const Box& tby = mfi.tilebox(By->ixType().toIntVect());
            const Box& tbz = mfi.tilebox(Bz->ixType().toIntVect());

            amrex::Array4<Real> const& sol_arr = m_solution[ilev]->array(mfi);
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

void run_ProjectionDivCleaner() {
    // Build Object, run, then delete
    ProjectionDivCleaner dc;
    dc.setSourceFromBfield();
    dc.solve();
    dc.correctBfield();
}
