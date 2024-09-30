/* Copyright 2023-2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *          S. Eric Clark (Helion Energy)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "HybridPICModel.H"

#include "EmbeddedBoundary/Enabled.H"
#include "Fields.H"
#include "Particles/MultiParticleContainer.H"
#include "WarpX.H"

using namespace amrex;
using warpx::fields::FieldType;

HybridPICModel::HybridPICModel ()
{
    ReadParameters();
}

void HybridPICModel::ReadParameters ()
{
    const ParmParse pp_hybrid("hybrid_pic_model");

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
    const bool n0_ref_given = utils::parser::queryWithParser(pp_hybrid, "n0_ref", m_n0_ref);
    if (m_gamma != 1.0 && !n0_ref_given) {
        Abort("hybrid_pic_model.n0_ref should be specified if hybrid_pic_model.gamma != 1");
    }

    pp_hybrid.query("plasma_resistivity(rho,J)", m_eta_expression);
    utils::parser::queryWithParser(pp_hybrid, "n_floor", m_n_floor);

    utils::parser::queryWithParser(pp_hybrid, "plasma_hyper_resistivity", m_eta_h);

    // convert electron temperature from eV to J
    m_elec_temp *= PhysConst::q_e;

    // external currents
    pp_hybrid.query("Jx_external_grid_function(x,y,z,t)", m_Jx_ext_grid_function);
    pp_hybrid.query("Jy_external_grid_function(x,y,z,t)", m_Jy_ext_grid_function);
    pp_hybrid.query("Jz_external_grid_function(x,y,z,t)", m_Jz_ext_grid_function);

    // external fields
    pp_hybrid.query("add_external_fields", m_add_external_fields);

    if (m_add_external_fields) {
        pp_hybrid.query("Bx_external_grid_function(x,y,z,t)", m_Bx_ext_grid_function);
        pp_hybrid.query("By_external_grid_function(x,y,z,t)", m_By_ext_grid_function);
        pp_hybrid.query("Bz_external_grid_function(x,y,z,t)", m_Bz_ext_grid_function);
        pp_hybrid.query("Ex_external_grid_function(x,y,z,t)", m_Ex_ext_grid_function);
        pp_hybrid.query("Ey_external_grid_function(x,y,z,t)", m_Ey_ext_grid_function);
        pp_hybrid.query("Ez_external_grid_function(x,y,z,t)", m_Ez_ext_grid_function);
    }
}

void HybridPICModel::AllocateLevelMFs (
    ablastr::fields::MultiFabRegister & fields,
    int lev, const BoxArray& ba, const DistributionMapping& dm,
    const int ncomps,
    const IntVect& ngJ, const IntVect& ngRho,
    const IntVect& ngEB,
    const IntVect& jx_nodal_flag,
    const IntVect& jy_nodal_flag,
    const IntVect& jz_nodal_flag,
    const IntVect& rho_nodal_flag,
    const IntVect& Ex_nodal_flag,
    const IntVect& Ey_nodal_flag,
    const IntVect& Ez_nodal_flag,
    const IntVect& Bx_nodal_flag,
    const IntVect& By_nodal_flag,
    const IntVect& Bz_nodal_flag)
{
    using ablastr::fields::Direction;

    // The "hybrid_electron_pressure_fp" multifab stores the electron pressure calculated
    // from the specified equation of state.
    // The "hybrid_rho_fp_temp" multifab is used to store the ion charge density
    // interpolated or extrapolated to appropriate timesteps.
    // The "hybrid_current_fp_temp" multifab is used to store the ion current density
    // interpolated or extrapolated to appropriate timesteps.
    // The "hybrid_current_fp_ampere" multifab stores the total current calculated as
    // the curl of B.
    fields.alloc_init(FieldType::hybrid_electron_pressure_fp,
        lev, amrex::convert(ba, rho_nodal_flag),
        dm, ncomps, ngRho, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_rho_fp_temp,
        lev, amrex::convert(ba, rho_nodal_flag),
        dm, ncomps, ngRho, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_temp, Direction{0},
        lev, amrex::convert(ba, jx_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_temp, Direction{1},
        lev, amrex::convert(ba, jy_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_temp, Direction{2},
        lev, amrex::convert(ba, jz_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);

    fields.alloc_init(FieldType::hybrid_current_fp_ampere, Direction{0},
        lev, amrex::convert(ba, jx_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_ampere, Direction{1},
        lev, amrex::convert(ba, jy_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_ampere, Direction{2},
        lev, amrex::convert(ba, jz_nodal_flag),
        dm, ncomps, ngJ, 0.0_rt);

    // the external current density multifab is made nodal to avoid needing to interpolate
    // to a nodal grid as has to be done for the ion and total current density multifabs
    // this also allows the external current multifab to not have any ghost cells
    fields.alloc_init(FieldType::hybrid_current_fp_external, Direction{0},
        lev, amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, IntVect(AMREX_D_DECL(0,0,0)), 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_external, Direction{1},
        lev, amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, IntVect(AMREX_D_DECL(0,0,0)), 0.0_rt);
    fields.alloc_init(FieldType::hybrid_current_fp_external, Direction{2},
        lev, amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, IntVect(AMREX_D_DECL(0,0,0)), 0.0_rt);

    if (m_add_external_fields) {
        // These are nodal to match when B-field is added in evaluation of Ohm's law
        fields.alloc_init(FieldType::hybrid_B_fp_external, Direction{0},
            lev, amrex::convert(ba, Bx_nodal_flag),
            dm, ncomps, ngEB, 0.0_rt);
        fields.alloc_init(FieldType::hybrid_B_fp_external, Direction{1},
            lev, amrex::convert(ba, By_nodal_flag),
            dm, ncomps, ngEB, 0.0_rt);
        fields.alloc_init(FieldType::hybrid_B_fp_external, Direction{2},
            lev, amrex::convert(ba, Bz_nodal_flag),
            dm, ncomps, ngEB, 0.0_rt);
        fields.alloc_init(FieldType::hybrid_E_fp_external, Direction{0},
            lev, amrex::convert(ba, Ex_nodal_flag),
            dm, ncomps, ngEB, 0.0_rt);
        fields.alloc_init(FieldType::hybrid_E_fp_external, Direction{1},
            lev, amrex::convert(ba, Ey_nodal_flag),
            dm, ncomps, ngEB, 0.0_rt);
        fields.alloc_init(FieldType::hybrid_E_fp_external, Direction{2},
            lev, amrex::convert(ba, Ez_nodal_flag),
            dm, ncomps, ngEB, 0.0_rt);
    }

#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (ncomps == 1),
        "Ohm's law solver only support m = 0 azimuthal mode at present.");
#endif
}

void HybridPICModel::InitData ()
{
    m_resistivity_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_eta_expression, {"rho","J"}));
    m_eta = m_resistivity_parser->compile<2>();
    const std::set<std::string> resistivity_symbols = m_resistivity_parser->symbols();
    m_resistivity_has_J_dependence += resistivity_symbols.count("J");

    m_J_external_parser[0] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Jx_ext_grid_function,{"x","y","z","t"}));
    m_J_external_parser[1] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Jy_ext_grid_function,{"x","y","z","t"}));
    m_J_external_parser[2] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Jz_ext_grid_function,{"x","y","z","t"}));
    m_J_external[0] = m_J_external_parser[0]->compile<4>();
    m_J_external[1] = m_J_external_parser[1]->compile<4>();
    m_J_external[2] = m_J_external_parser[2]->compile<4>();

    // check if the external current parsers depend on time
    for (int i=0; i<3; i++) {
        const std::set<std::string> J_ext_symbols = m_J_external_parser[i]->symbols();
        m_external_field_has_time_dependence += J_ext_symbols.count("t");
    }

    auto & warpx = WarpX::GetInstance();

    m_B_external_parser[0] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Bx_ext_grid_function,{"x","y","z","t"}));
    m_B_external_parser[1] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_By_ext_grid_function,{"x","y","z","t"}));
    m_B_external_parser[2] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Bz_ext_grid_function,{"x","y","z","t"}));
    m_B_external[0] = m_B_external_parser[0]->compile<4>();
    m_B_external[1] = m_B_external_parser[1]->compile<4>();
    m_B_external[2] = m_B_external_parser[2]->compile<4>();

    m_E_external_parser[0] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Ex_ext_grid_function,{"x","y","z","t"}));
    m_E_external_parser[1] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Ey_ext_grid_function,{"x","y","z","t"}));
    m_E_external_parser[2] = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_Ez_ext_grid_function,{"x","y","z","t"}));
    m_E_external[0] = m_E_external_parser[0]->compile<4>();
    m_E_external[1] = m_E_external_parser[1]->compile<4>();
    m_E_external[2] = m_E_external_parser[2]->compile<4>();
    using ablastr::fields::Direction;

    // Get the grid staggering of the fields involved in calculating E
    amrex::IntVect Jx_stag = warpx.m_fields.get(FieldType::current_fp, Direction{0}, 0)->ixType().toIntVect();
    amrex::IntVect Jy_stag = warpx.m_fields.get(FieldType::current_fp, Direction{1}, 0)->ixType().toIntVect();
    amrex::IntVect Jz_stag = warpx.m_fields.get(FieldType::current_fp, Direction{2}, 0)->ixType().toIntVect();
    amrex::IntVect Bx_stag = warpx.m_fields.get(FieldType::Bfield_fp, Direction{0}, 0)->ixType().toIntVect();
    amrex::IntVect By_stag = warpx.m_fields.get(FieldType::Bfield_fp, Direction{1}, 0)->ixType().toIntVect();
    amrex::IntVect Bz_stag = warpx.m_fields.get(FieldType::Bfield_fp, Direction{2}, 0)->ixType().toIntVect();
    amrex::IntVect Ex_stag = warpx.m_fields.get(FieldType::Efield_fp, Direction{0}, 0)->ixType().toIntVect();
    amrex::IntVect Ey_stag = warpx.m_fields.get(FieldType::Efield_fp, Direction{1}, 0)->ixType().toIntVect();
    amrex::IntVect Ez_stag = warpx.m_fields.get(FieldType::Efield_fp, Direction{2}, 0)->ixType().toIntVect();

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

    // Initialize external current - note that this approach skips the check
    // if the current is time dependent which is what needs to be done to
    // write time independent fields on the first step.
    GetCurrentExternal(true);
    if (m_add_external_fields)
        GetFieldsExternal(warpx.gett_new(0));
}

void HybridPICModel::GetCurrentExternal (bool skip_check /*=false*/)
{
    if (!skip_check && !m_external_field_has_time_dependence) { return; }

    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        GetExternalFieldFromExpression(FieldType::hybrid_current_fp_external, m_J_external, lev);
    }
}

void HybridPICModel::GetFieldsExternal (amrex::Real t)
{
    if (!m_add_external_fields) return;

    using ablastr::fields::Direction;
    auto& warpx = WarpX::GetInstance();

    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        GetExternalFieldFromExpression(
            FieldType::hybrid_B_fp_external,
            m_B_external, lev, t);
        GetExternalFieldFromExpression(
            FieldType::hybrid_E_fp_external,
            m_E_external, lev, t);
        for (int idim=0; idim < 3; idim++) {
            auto mf_Bext = warpx.m_fields.get(FieldType::hybrid_B_fp_external, Direction{idim}, lev);
            mf_Bext->FillBoundary(warpx.Geom(lev).periodicity());
            auto mf_Eext = warpx.m_fields.get(FieldType::hybrid_E_fp_external, Direction{idim}, lev);
            mf_Eext->FillBoundary(warpx.Geom(lev).periodicity());
        }
    }
}

void HybridPICModel::GetExternalFieldFromExpression (
    FieldType field_type,
    std::array< amrex::ParserExecutor<4>, 3> const& expression,
    int lev)
{
    auto & warpx = WarpX::GetInstance();
    GetExternalFieldFromExpression(field_type, expression, lev, warpx.gett_new(lev));
}

void HybridPICModel::GetExternalFieldFromExpression (
    FieldType field_type,
    std::array< amrex::ParserExecutor<4>, 3> const& expression,
    int lev, amrex::Real t)
{
    using ablastr::fields::Direction;

    // This logic matches closely to WarpX::InitializeExternalFieldsOnGridUsingParser
    // except that the parsers include time dependence.
    auto & warpx = WarpX::GetInstance();
    auto dx_lev = warpx.Geom(lev).CellSizeArray();
    const RealBox& real_box = warpx.Geom(lev).ProbDomain();

    ablastr::fields::VectorField field = warpx.m_fields.get_alldirs(field_type, lev);

    amrex::MultiFab* mfx = field[0];
    amrex::MultiFab* mfy = field[1];
    amrex::MultiFab* mfz = field[2];

    const amrex::IntVect x_nodal_flag = mfx->ixType().toIntVect();
    const amrex::IntVect y_nodal_flag = mfy->ixType().toIntVect();
    const amrex::IntVect z_nodal_flag = mfz->ixType().toIntVect();

    // avoid implicit lambda capture
    auto x_external = expression[0];
    auto y_external = expression[1];
    auto z_external = expression[2];

    for ( MFIter mfi(*mfx, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& tbx = mfi.tilebox( x_nodal_flag, mfx->nGrowVect() );
        const amrex::Box& tby = mfi.tilebox( y_nodal_flag, mfy->nGrowVect() );
        const amrex::Box& tbz = mfi.tilebox( z_nodal_flag, mfz->nGrowVect() );

        auto const& mfxfab = mfx->array(mfi);
        auto const& mfyfab = mfy->array(mfi);
        auto const& mfzfab = mfz->array(mfi);

        amrex::Box lxb, lyb, lzb;
        amrex::Array4<const amrex::Real> lx, ly, lz;
        if (EB::enabled()) {

            auto const& mf_lx = *warpx.m_fields.get(FieldType::edge_lengths, Direction{0}, lev);
            auto const& mf_ly = *warpx.m_fields.get(FieldType::edge_lengths, Direction{1}, lev);
            auto const& mf_lz = *warpx.m_fields.get(FieldType::edge_lengths, Direction{2}, lev);

            lxb = mfi.growntilebox(mf_lx.nGrowVect());
            lyb = mfi.growntilebox(mf_ly.nGrowVect());
            lzb = mfi.growntilebox(mf_lz.nGrowVect());

            lx = mf_lx.array(mfi);
            ly = mf_ly.array(mfi);
            lz = mf_lz.array(mfi);
        }

        amrex::ParallelFor (tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // skip if node is covered by an embedded boundary or outside of lx box array
                if (lx && (!lxb.contains({i, j, k}) || lx(i, j, k) <= 0)) { return; }

                // Shift required in the x-, y-, or z- position
                // depending on the index type of the multifab
#if defined(WARPX_DIM_1D_Z)
                const amrex::Real x = 0._rt;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - x_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real z = j*dx_lev[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const amrex::Real fac_x = (1._rt - x_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - x_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#else
                const amrex::Real fac_x = (1._rt - x_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real fac_y = (1._rt - x_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                const amrex::Real fac_z = (1._rt - x_nodal_flag[2]) * dx_lev[2] * 0.5_rt;
                const amrex::Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the x-component of the field.
                mfxfab(i,j,k) = x_external(x,y,z,t);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // skip if node is covered by an embedded boundary or outside of ly box array
                if (ly && (!lyb.contains({i, j, k}) || ly(i, j, k) <= 0)) { return; }

#if defined(WARPX_DIM_1D_Z)
                const amrex::Real x = 0._rt;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - y_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real z = j*dx_lev[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const amrex::Real fac_x = (1._rt - y_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - y_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif defined(WARPX_DIM_3D)
                const amrex::Real fac_x = (1._rt - y_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real fac_y = (1._rt - y_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                const amrex::Real fac_z = (1._rt - y_nodal_flag[2]) * dx_lev[2] * 0.5_rt;
                const amrex::Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the y-component of the field.
                mfyfab(i,j,k)  = y_external(x,y,z,t);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // skip if node is covered by an embedded boundary or outside of lz box array
                if (lz && (!lzb.contains({i, j, k}) || lz(i, j, k) <= 0)) { return; }

#if defined(WARPX_DIM_1D_Z)
                const amrex::Real x = 0._rt;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - z_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real z = j*dx_lev[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const amrex::Real fac_x = (1._rt - z_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real y = 0._rt;
                const amrex::Real fac_z = (1._rt - z_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif defined(WARPX_DIM_3D)
                const amrex::Real fac_x = (1._rt - z_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                const amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                const amrex::Real fac_y = (1._rt - z_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                const amrex::Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                const amrex::Real fac_z = (1._rt - z_nodal_flag[2]) * dx_lev[2] * 0.5_rt;
                const amrex::Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the z-component of the field.
                mfzfab(i,j,k) = z_external(x,y,z,t);
            }
        );
    }
    amrex::Gpu::streamSynchronize();
}

void HybridPICModel::CalculateCurrentAmpere (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths)
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        CalculateCurrentAmpere(Bfield[lev], edge_lengths[lev], lev);
    }
}

void HybridPICModel::CalculateCurrentAmpere (
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& edge_lengths,
    const int lev)
{
    WARPX_PROFILE("WarpX::CalculateCurrentAmpere()");

    auto& warpx = WarpX::GetInstance();
    ablastr::fields::VectorField current_fp_ampere = warpx.m_fields.get_alldirs(FieldType::hybrid_current_fp_ampere, lev);
    warpx.get_pointer_fdtd_solver_fp(lev)->CalculateCurrentAmpere(
        current_fp_ampere, Bfield, edge_lengths, lev
    );

    // we shouldn't apply the boundary condition to J since J = J_i - J_e but
    // the boundary correction was already applied to J_i and the B-field
    // boundary ensures that J itself complies with the boundary conditions, right?
    // ApplyJfieldBoundary(lev, Jfield[0].get(), Jfield[1].get(), Jfield[2].get());
    for (int i=0; i<3; i++) { current_fp_ampere[i]->FillBoundary(warpx.Geom(lev).periodicity()); }
}

void HybridPICModel::HybridPICSolveE (
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths,
    amrex::Real t,
    const bool solve_for_Faraday) const
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        HybridPICSolveE(
            Efield[lev], Jfield[lev], Bfield[lev], *rhofield[lev],
            edge_lengths[lev], t, lev, solve_for_Faraday
        );
    }
}

void HybridPICModel::HybridPICSolveE (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Jfield,
    ablastr::fields::VectorField const& Bfield,
    amrex::MultiFab const& rhofield,
    ablastr::fields::VectorField const& edge_lengths,
    amrex::Real t,
    const int lev, const bool solve_for_Faraday) const
{
    WARPX_PROFILE("WarpX::HybridPICSolveE()");

    HybridPICSolveE(
        Efield, Jfield, Bfield, rhofield, edge_lengths, t, lev,
        PatchType::fine, solve_for_Faraday
    );
    if (lev > 0)
    {
        amrex::Abort(Utils::TextMsg::Err(
        "HybridPICSolveE: Only one level implemented for hybrid-PIC solver."));
    }
}

void HybridPICModel::HybridPICSolveE (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Jfield,
    ablastr::fields::VectorField const& Bfield,
    amrex::MultiFab const& rhofield,
    ablastr::fields::VectorField const& edge_lengths,
    amrex::Real t,
    const int lev, PatchType patch_type,
    const bool solve_for_Faraday) const
{

    auto& warpx = WarpX::GetInstance();

    ablastr::fields::VectorField current_fp_ampere = warpx.m_fields.get_alldirs(FieldType::hybrid_current_fp_ampere, lev);
    const ablastr::fields::VectorField current_fp_external = warpx.m_fields.get_alldirs(FieldType::hybrid_current_fp_external, lev);
    const ablastr::fields::ScalarField electron_pressure_fp = warpx.m_fields.get(FieldType::hybrid_electron_pressure_fp, lev);

    // Solve E field in regular cells
    warpx.get_pointer_fdtd_solver_fp(lev)->HybridPICSolveE(
        Efield, current_fp_ampere, Jfield, current_fp_external,
        Bfield, rhofield,
        *electron_pressure_fp,
        edge_lengths, t, lev, this, solve_for_Faraday
    );
    warpx.ApplyEfieldBoundary(lev, patch_type);
}

void HybridPICModel::CalculateElectronPressure() const
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        CalculateElectronPressure(lev);
    }
}

void HybridPICModel::CalculateElectronPressure(const int lev) const
{
    WARPX_PROFILE("WarpX::CalculateElectronPressure()");

    auto& warpx = WarpX::GetInstance();
    ablastr::fields::ScalarField electron_pressure_fp = warpx.m_fields.get(FieldType::hybrid_electron_pressure_fp, lev);
    ablastr::fields::ScalarField rho_fp = warpx.m_fields.get(FieldType::rho_fp, lev);

    // Calculate the electron pressure using rho^{n+1}.
    FillElectronPressureMF(
        *electron_pressure_fp,
        *rho_fp
    );
    warpx.ApplyElectronPressureBoundary(lev, PatchType::fine);
    electron_pressure_fp->FillBoundary(warpx.Geom(lev).periodicity());
}

void HybridPICModel::FillElectronPressureMF (
    amrex::MultiFab& Pe_field,
    amrex::MultiFab const& rho_field
) const
{
    const auto n0_ref = m_n0_ref;
    const auto elec_temp = m_elec_temp;
    const auto gamma = m_gamma;

    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(Pe_field, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Extract field data for this grid/tile
        Array4<Real const> const& rho = rho_field.const_array(mfi);
        Array4<Real> const& Pe = Pe_field.array(mfi);

        // Extract tileboxes for which to loop
        const Box& tilebox  = mfi.tilebox();

        ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Pe(i, j, k) = ElectronPressure::get_pressure(
                n0_ref, elec_temp, gamma, rho(i, j, k)
            );
        });
    }
}

void HybridPICModel::BfieldEvolveRK (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField  const& edge_lengths,
    amrex::Real t, amrex::Real dt, DtType dt_type,
    IntVect ng, std::optional<bool> nodal_sync )
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev)
    {
        BfieldEvolveRK(
            Bfield, Efield, Jfield, rhofield, edge_lengths, t, dt, lev, dt_type,
            ng, nodal_sync
        );
    }
}

void HybridPICModel::BfieldEvolveRK (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths,
    amrex::Real t, amrex::Real dt, int lev, DtType dt_type,
    IntVect ng, std::optional<bool> nodal_sync )
{
    // Make copies of the B-field multifabs at t = n and create multifabs for
    // each direction to store the Runge-Kutta intermediate terms. Each
    // multifab has 2 components for the different terms that need to be stored.
    std::array< MultiFab, 3 > B_old;
    std::array< MultiFab, 3 > K;
    for (int ii = 0; ii < 3; ii++)
    {
        B_old[ii] = MultiFab(
            Bfield[lev][ii]->boxArray(), Bfield[lev][ii]->DistributionMap(), 1,
            Bfield[lev][ii]->nGrowVect()
        );
        MultiFab::Copy(B_old[ii], *Bfield[lev][ii], 0, 0, 1, ng);

        K[ii] = MultiFab(
            Bfield[lev][ii]->boxArray(), Bfield[lev][ii]->DistributionMap(), 2,
            Bfield[lev][ii]->nGrowVect()
        );
        K[ii].setVal(0.0);
    }

    amrex::Real t_eval = t;

    // The Runge-Kutta scheme begins here.
    // Step 1:
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        t_eval, 0.5_rt*dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + 0.5 * dt * [-curl x E(B_old)] = B_old + 0.5 * dt * K0.
    for (int ii = 0; ii < 3; ii++)
    {
        // Extract 0.5 * dt * K0 for each direction into index 0 of K.
        MultiFab::LinComb(
            K[ii], 1._rt, *Bfield[lev][ii], 0, -1._rt, B_old[ii], 0, 0, 1, ng
        );
    }

    // Step 2:
    t_eval = t+0.5_rt*dt;
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        t_eval, 0.5_rt*dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + 0.5 * dt * K0 + 0.5 * dt * [-curl x E(B_old + 0.5 * dt * K1)]
    //       = B_old + 0.5 * dt * K0 + 0.5 * dt * K1
    for (int ii = 0; ii < 3; ii++)
    {
        // Subtract 0.5 * dt * K0 from the Bfield for each direction, to get
        // B_new = B_old + 0.5 * dt * K1.
        MultiFab::Subtract(*Bfield[lev][ii], K[ii], 0, 0, 1, ng);
        // Extract 0.5 * dt * K1 for each direction into index 1 of K.
        MultiFab::LinComb(
            K[ii], 1._rt, *Bfield[lev][ii], 0, -1._rt, B_old[ii], 0, 1, 1, ng
        );
    }

    // Step 3:
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        t_eval, dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + 0.5 * dt * K1 + dt * [-curl  x E(B_old + 0.5 * dt * K1)]
    //       = B_old + 0.5 * dt * K1 + dt * K2
    for (int ii = 0; ii < 3; ii++)
    {
        // Subtract 0.5 * dt * K1 from the Bfield for each direction to get
        // B_new = B_old + dt * K2.
        MultiFab::Subtract(*Bfield[lev][ii], K[ii], 1, 0, 1, ng);
    }

    // Step 4:
    t_eval = t + dt;
    FieldPush(
        Bfield, Efield, Jfield, rhofield, edge_lengths,
        t_eval, 0.5_rt*dt, dt_type, ng, nodal_sync
    );

    // The Bfield is now given by:
    // B_new = B_old + dt * K2 + 0.5 * dt * [-curl x E(B_old + dt * K2)]
    //       = B_old + dt * K2 + 0.5 * dt * K3
    for (int ii = 0; ii < 3; ii++)
    {
        // Subtract B_old from the Bfield for each direction, to get
        // B = dt * K2 + 0.5 * dt * K3.
        MultiFab::Subtract(*Bfield[lev][ii], B_old[ii], 0, 0, 1, ng);

        // Add dt * K2 + 0.5 * dt * K3 to index 0 of K (= 0.5 * dt * K0).
        MultiFab::Add(K[ii], *Bfield[lev][ii], 0, 0, 1, ng);

        // Add 2 * 0.5 * dt * K1 to index 0 of K.
        MultiFab::LinComb(
            K[ii], 1.0, K[ii], 0, 2.0, K[ii], 1, 0, 1, ng
        );

        // Overwrite the Bfield with the Runge-Kutta sum:
        // B_new = B_old + 1/3 * dt * (0.5 * K0 + K1 + K2 + 0.5 * K3).
        MultiFab::LinComb(
            *Bfield[lev][ii], 1.0, B_old[ii], 0, 1.0/3.0, K[ii], 0, 0, 1, ng
        );
    }
}

void HybridPICModel::FieldPush (
    ablastr::fields::MultiLevelVectorField const& Bfield,
    ablastr::fields::MultiLevelVectorField const& Efield,
    ablastr::fields::MultiLevelVectorField const& Jfield,
    ablastr::fields::MultiLevelScalarField const& rhofield,
    ablastr::fields::MultiLevelVectorField const& edge_lengths,
    amrex::Real t, amrex::Real dt, DtType dt_type,
    IntVect ng, std::optional<bool> nodal_sync )
{
    auto& warpx = WarpX::GetInstance();

    // Calculate J = curl x B / mu0
    CalculateCurrentAmpere(Bfield, edge_lengths);
    // Calculate the E-field from Ohm's law
    HybridPICSolveE(Efield, Jfield, Bfield, rhofield, edge_lengths, t, true);
    warpx.FillBoundaryE(ng, nodal_sync);
    // Push forward the B-field using Faraday's law
    warpx.EvolveB(dt, dt_type);
    warpx.FillBoundaryB(ng, nodal_sync);
}
