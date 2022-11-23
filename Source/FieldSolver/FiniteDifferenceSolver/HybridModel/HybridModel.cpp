/* Copyright 2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald
 *
 * License: BSD-3-Clause-LBNL
 */

#include "HybridModel.H"

#include "Utils/WarpXConst.H"
#include "Utils/Parser/ParserUtils.H"
#include "WarpX.H"


using namespace amrex;

HybridModel::HybridModel ()
{
    ReadParameters();
}

void
HybridModel::ReadParameters ()
{
    ParmParse pp_hybrid("hybridmodel");

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
    amrex::IntVect rho_stag = warpx.getrho_fp(0).ixType().toIntVect();
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
        rho_IndexType[idim]   = rho_stag[idim];
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
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_1D_Z)
        rho_IndexType[2]   = 0;
        Jx_IndexType[2]    = 0;
        Jy_IndexType[2]    = 0;
        Jz_IndexType[2]    = 0;
        Bx_IndexType[2]    = 0;
        By_IndexType[2]    = 0;
        Bz_IndexType[2]    = 0;
        Ex_IndexType[2]    = 0;
        Ey_IndexType[2]    = 0;
        Ez_IndexType[2]    = 0;
#endif
#if defined(WARPX_DIM_1D_Z)
        rho_IndexType[1]   = 0;
        Jx_IndexType[1]    = 0;
        Jy_IndexType[1]    = 0;
        Jz_IndexType[1]    = 0;
        Bx_IndexType[1]    = 0;
        By_IndexType[1]    = 0;
        Bz_IndexType[1]    = 0;
        Ex_IndexType[1]    = 0;
        Ey_IndexType[1]    = 0;
        Ez_IndexType[1]    = 0;
#endif
}