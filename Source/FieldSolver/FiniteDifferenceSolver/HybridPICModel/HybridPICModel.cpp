/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "HybridPICModel.H"

using namespace amrex;

HybridPICModel::HybridPICModel ( int nlevs_max )
{
    ReadParameters();
    AllocateMFs(nlevs_max);
}

void HybridPICModel::ReadParameters ()
{
    ParmParse pp_hybrid("hybrid_pic_model");

    // The B-field update is subcycled to improve stability - the number
    // of sub steps can be specified by the user (defaults to 50).
    utils::parser::queryWithParser(pp_hybrid, "substeps", m_substeps);

    // The hybrid model requires an electron temperature, reference density
    // and exponent to be given. These values will be used to calculate the
    // electron pressure according to p = n0 * Te * (n/n0)^gamma
    utils::parser::queryWithParser(pp_hybrid, "gamma", m_gamma);
    if (!utils::parser::queryWithParser(pp_hybrid, "elec_temp", m_elec_temp)) {
        Abort("hybrid_pic_model.elec_temp must be specified when using the hybrid solver");
    }
    bool n0_ref_given = utils::parser::queryWithParser(pp_hybrid, "n0_ref", m_n0_ref);
    if (m_gamma != 1.0 && !n0_ref_given) {
        Abort("hybrid_pic_model.n0_ref should be specified if hybrid_pic_model.gamma != 1");
    }

    pp_hybrid.query("plasma_resistivity(rho)", m_eta_expression);
    utils::parser::queryWithParser(pp_hybrid, "n_floor", m_n_floor);

    // convert electron temperature from eV to J
    m_elec_temp *= PhysConst::q_e;
}

void HybridPICModel::AllocateMFs (int nlevs_max)
{
    electron_pressure_fp.resize(nlevs_max);
    rho_fp_temp.resize(nlevs_max);
    current_fp_temp.resize(nlevs_max);
    current_fp_ampere.resize(nlevs_max);
}

void HybridPICModel::AllocateLevelMFs (int lev, const BoxArray& ba, const DistributionMapping& dm,
                                       const int ncomps, const IntVect& ngJ, const IntVect& ngRho,
                                       const IntVect& jx_nodal_flag,
                                       const IntVect& jy_nodal_flag,
                                       const IntVect& jz_nodal_flag,
                                       const IntVect& rho_nodal_flag)
{
    // The "electron_pressure_fp" multifab stores the electron pressure calculated
    // from the specified equation of state.
    // The "rho_fp_temp" multifab is used to store the ion charge density
    // interpolated or extrapolated to appropriate timesteps.
    // The "current_fp_temp" multifab is used to store the ion current density
    // interpolated or extrapolated to appropriate timesteps.
    // The "current_fp_ampere" multifab stores the total current calculated as
    // the curl of B.
    WarpX::AllocInitMultiFab(electron_pressure_fp[lev], amrex::convert(ba, rho_nodal_flag),
        dm, ncomps, ngRho, lev, "electron_pressure_fp", 0.0_rt);

    WarpX::AllocInitMultiFab(rho_fp_temp[lev], amrex::convert(ba, rho_nodal_flag),
        dm, ncomps, ngRho, lev, "rho_fp_temp", 0.0_rt);

    WarpX::AllocInitMultiFab(current_fp_temp[lev][0], amrex::convert(ba, jx_nodal_flag),
        dm, ncomps, ngJ, lev, "current_fp_temp[x]", 0.0_rt);
    WarpX::AllocInitMultiFab(current_fp_temp[lev][1], amrex::convert(ba, jy_nodal_flag),
        dm, ncomps, ngJ, lev, "current_fp_temp[y]", 0.0_rt);
    WarpX::AllocInitMultiFab(current_fp_temp[lev][2], amrex::convert(ba, jz_nodal_flag),
        dm, ncomps, ngJ, lev, "current_fp_temp[z]", 0.0_rt);

    WarpX::AllocInitMultiFab(current_fp_ampere[lev][0], amrex::convert(ba, jx_nodal_flag),
        dm, ncomps, ngJ, lev, "current_fp_ampere[x]", 0.0_rt);
    WarpX::AllocInitMultiFab(current_fp_ampere[lev][1], amrex::convert(ba, jy_nodal_flag),
        dm, ncomps, ngJ, lev, "current_fp_ampere[y]", 0.0_rt);
    WarpX::AllocInitMultiFab(current_fp_ampere[lev][2], amrex::convert(ba, jz_nodal_flag),
        dm, ncomps, ngJ, lev, "current_fp_ampere[z]", 0.0_rt);
}

void HybridPICModel::ClearLevel (int lev)
{
    electron_pressure_fp[lev].reset();
    rho_fp_temp[lev].reset();
    for (int i = 0; i < 3; ++i) {
        current_fp_temp[lev][i].reset();
        current_fp_ampere[lev][i].reset();
    }
}

void HybridPICModel::InitData ()
{
    m_resistivity_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_eta_expression, {"rho"}));
    m_eta = m_resistivity_parser->compile<1>();

    auto & warpx = WarpX::GetInstance();

    // Get the grid staggering of the fields involved in calculating E
    amrex::IntVect Jx_stag = warpx.getcurrent_fp(0,0).ixType().toIntVect();
    amrex::IntVect Jy_stag = warpx.getcurrent_fp(0,1).ixType().toIntVect();
    amrex::IntVect Jz_stag = warpx.getcurrent_fp(0,2).ixType().toIntVect();
    amrex::IntVect Bx_stag = warpx.getBfield_fp(0,0).ixType().toIntVect();
    amrex::IntVect By_stag = warpx.getBfield_fp(0,1).ixType().toIntVect();
    amrex::IntVect Bz_stag = warpx.getBfield_fp(0,2).ixType().toIntVect();
    amrex::IntVect Ex_stag = warpx.getEfield_fp(0,0).ixType().toIntVect();
    amrex::IntVect Ey_stag = warpx.getEfield_fp(0,1).ixType().toIntVect();
    amrex::IntVect Ez_stag = warpx.getEfield_fp(0,2).ixType().toIntVect();

    // Check that the grid types are appropriate
    const bool appropriate_grids = (
#if   defined(WARPX_DIM_1D_Z)
        // AMReX convention: x = missing dimension, y = missing dimension, z = only dimension
        Ex_stag == IntVect(1) && Ey_stag == IntVect(1) && Ez_stag == IntVect(0) &&
        Bx_stag == IntVect(0) && By_stag == IntVect(0) && Bz_stag == IntVect(1) &&
#elif   defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        // AMReX convention: x = first dimension, y = missing dimension, z = second dimension
        Ex_stag == IntVect(0,1) && Ey_stag == IntVect(1,1) && Ez_stag == IntVect(1,0) &&
        Bx_stag == IntVect(1,0) && By_stag == IntVect(0,0) && Bz_stag == IntVect(0,1) &&
#elif defined(WARPX_DIM_3D)
        Ex_stag == IntVect(0,1,1) && Ey_stag == IntVect(1,0,1) && Ez_stag == IntVect(1,1,0) &&
        Bx_stag == IntVect(1,0,0) && By_stag == IntVect(0,1,0) && Bz_stag == IntVect(0,0,1) &&
#endif
        Jx_stag == Ex_stag && Jy_stag == Ey_stag && Jz_stag == Ez_stag
    );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        appropriate_grids,
        "Ohm's law E-solve only works with staggered (Yee) grids.");

    // copy data to device
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Jx_IndexType[idim]    = Jx_stag[idim];
        Jy_IndexType[idim]    = Jy_stag[idim];
        Jz_IndexType[idim]    = Jz_stag[idim];
        Bx_IndexType[idim]    = Bx_stag[idim];
        By_IndexType[idim]    = By_stag[idim];
        Bz_IndexType[idim]    = Bz_stag[idim];
        Ex_IndexType[idim]    = Ex_stag[idim];
        Ey_IndexType[idim]    = Ey_stag[idim];
        Ez_IndexType[idim]    = Ez_stag[idim];
    }

    // Below we set all the unused dimensions to have nodal values for J, B & E
    // since these values will be interpolated onto a nodal grid - if this is
    // not done the Interp function returns nonsense values.
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_1D_Z)
    Jx_IndexType[2]    = 1;
    Jy_IndexType[2]    = 1;
    Jz_IndexType[2]    = 1;
    Bx_IndexType[2]    = 1;
    By_IndexType[2]    = 1;
    Bz_IndexType[2]    = 1;
    Ex_IndexType[2]    = 1;
    Ey_IndexType[2]    = 1;
    Ez_IndexType[2]    = 1;
#endif
#if defined(WARPX_DIM_1D_Z)
    Jx_IndexType[1]    = 1;
    Jy_IndexType[1]    = 1;
    Jz_IndexType[1]    = 1;
    Bx_IndexType[1]    = 1;
    By_IndexType[1]    = 1;
    Bz_IndexType[1]    = 1;
    Ex_IndexType[1]    = 1;
    Ey_IndexType[1]    = 1;
    Ez_IndexType[1]    = 1;
#endif
}

void HybridPICModel::CalculateCurrentAmpere (
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Bfield,
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths)
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        CalculateCurrentAmpere(Bfield[lev], edge_lengths[lev], lev);
    }
}

void HybridPICModel::CalculateCurrentAmpere (
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
    const int lev)
{
    WARPX_PROFILE("WarpX::CalculateCurrentAmpere()");

    auto& warpx = WarpX::GetInstance();
    warpx.get_pointer_fdtd_solver_fp(lev)->CalculateCurrentAmpere(
        current_fp_ampere[lev], Bfield, edge_lengths, lev
    );

    // we shouldn't apply the boundary condition to J since J = J_i - J_e but
    // the boundary correction was already applied to J_i and the B-field
    // boundary ensures that J itself complies with the boundary conditions, right?
    // ApplyJfieldBoundary(lev, Jfield[0].get(), Jfield[1].get(), Jfield[2].get());
    for (int i=0; i<3; i++) current_fp_ampere[lev][i]->FillBoundary(warpx.Geom(lev).periodicity());
}

void HybridPICModel::HybridPICSolveE (
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> & Efield,
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Jfield,
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Bfield,
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& rhofield,
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths,
    const bool include_resistivity_term)
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        HybridPICSolveE(
            Efield[lev], Jfield[lev], Bfield[lev], rhofield[lev],
            edge_lengths[lev], lev, include_resistivity_term
        );
    }
}

void HybridPICModel::HybridPICSolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3> & Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
    const int lev, const bool include_resistivity_term)
{
    WARPX_PROFILE("WarpX::HybridPICSolveE()");

    HybridPICSolveE(
        Efield, Jfield, Bfield, rhofield, edge_lengths, lev,
        PatchType::fine, include_resistivity_term
    );
    if (lev > 0)
    {
        amrex::Abort(Utils::TextMsg::Err(
        "HybridPICSolveE: Only one level implemented for hybrid-PIC solver."));
    }
}

void HybridPICModel::HybridPICSolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3> & Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Jfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
    const int lev, PatchType patch_type,
    const bool include_resistivity_term)
{
    auto& warpx = WarpX::GetInstance();

    // Solve E field in regular cells
    warpx.get_pointer_fdtd_solver_fp(lev)->HybridPICSolveE(
        Efield, current_fp_ampere[lev], Jfield, Bfield, rhofield,
        electron_pressure_fp[lev],
        edge_lengths, lev, this, include_resistivity_term
    );
    warpx.ApplyEfieldBoundary(lev, patch_type);
}

void HybridPICModel::CalculateElectronPressure(DtType a_dt_type)
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        CalculateElectronPressure(lev, a_dt_type);
    }
}

void HybridPICModel::CalculateElectronPressure(const int lev, DtType a_dt_type)
{
    WARPX_PROFILE("WarpX::CalculateElectronPressure()");

    auto& warpx = WarpX::GetInstance();
    // The full step uses rho^{n+1}, otherwise use the old or averaged
    // charge density.
    if (a_dt_type == DtType::Full) {
        FillElectronPressureMF(
            electron_pressure_fp[lev], warpx.get_pointer_rho_fp(lev)
        );
    } else {
        FillElectronPressureMF(
            electron_pressure_fp[lev], rho_fp_temp[lev].get()
        );
    }
    warpx.ApplyElectronPressureBoundary(lev, PatchType::fine);
    electron_pressure_fp[lev]->FillBoundary(warpx.Geom(lev).periodicity());
}

void HybridPICModel::FillElectronPressureMF (
    std::unique_ptr<amrex::MultiFab> const& Pe_field,
    amrex::MultiFab* const& rho_field )
{
    const auto n0_ref = m_n0_ref;
    const auto elec_temp = m_elec_temp;
    const auto gamma = m_gamma;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Pe_field, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Extract field data for this grid/tile
        Array4<Real const> const& rho = rho_field->const_array(mfi);
        Array4<Real> const& Pe = Pe_field->array(mfi);

        // Extract tileboxes for which to loop
        const Box& tilebox  = mfi.tilebox();

        ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Pe(i, j, k) = ElectronPressure::get_pressure(
                n0_ref, elec_temp, gamma, rho(i, j, k)
            );
        });
    }
}
