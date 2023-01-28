/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "HybridModel.H"
#include "Utils/Parser/ParserUtils.H"


using namespace amrex;

HybridModel::HybridModel ()
{
    ReadParameters();
}

void
HybridModel::ReadParameters ()
{
    ParmParse pp_hybrid("hybridmodel");

    // The B-field update is subcycled to improve stability - the number
    // of sub steps can be specified by the user (defaults to 50).
    utils::parser::queryWithParser(pp_hybrid, "substeps", m_substeps);

    // The hybrid model requires an electron temperature, reference density
    // and exponent to be given. These values will be used to calculate the
    // electron pressure according to p = n0 * Te * (n/n0)^gamma
    if (!utils::parser::queryWithParser(pp_hybrid, "elec_temp", m_elec_temp)) {
        Abort("hybridmodel.elec_temp must be specified when using the hybrid solver");
    }
    if (!utils::parser::queryWithParser(pp_hybrid, "n0_ref", m_n0_ref)) {
        Abort("hybridmodel.n0_ref must be specified when using the hybrid solver");
    }
    utils::parser::queryWithParser(pp_hybrid, "gamma", m_gamma);
    utils::parser::queryWithParser(pp_hybrid, "plasma_resistivity", m_eta);

    // convert electron temperature from eV to J
    m_elec_temp *= PhysConst::q_e;
}

void
HybridModel::InitData ()
{
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

void HybridModel::FillElectronPressureMF (
    std::unique_ptr<amrex::MultiFab> const& Pe_field,
    std::unique_ptr<amrex::MultiFab> const& rho_field,
    DtType a_dt_type )
{
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

            Real rho_val;
            if (a_dt_type == DtType::FirstHalf) {
                // use rho^{n}
                rho_val = rho(i, j, k, 0);
            }
            else if (a_dt_type == DtType::SecondHalf) {
                // use rho^{n+1/2}
                rho_val = 0.5 * (rho(i, j, k, 0) + rho(i, j, k, 1));
            }
            else {
                // use rho^{n+1}
                rho_val = rho(i, j, k, 1);
            }

            Pe(i, j, k) = ElectronPressure::get_pressure(
                m_n0_ref, m_elec_temp, m_gamma, rho_val
            );
        });
    }
}