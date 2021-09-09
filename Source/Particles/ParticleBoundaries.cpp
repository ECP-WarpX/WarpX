/* Copyright 2021 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleBoundaries.H"

ParticleBoundaries::ParticleBoundaries () noexcept
{
    SetAll(ParticleBoundaryType::Absorbing);
    reflect_all_velocities = false;
}

void
ParticleBoundaries::Set_reflect_all_velocities (bool flag)
{
    reflect_all_velocities = flag;
}

void
ParticleBoundaries::SetAll (ParticleBoundaryType bc)
{
    xmin_bc = bc;
    xmax_bc = bc;
    ymin_bc = bc;
    ymax_bc = bc;
    zmin_bc = bc;
    zmax_bc = bc;
}

void
ParticleBoundaries::SetBoundsX (ParticleBoundaryType bc_lo, ParticleBoundaryType bc_hi)
{
    xmin_bc = bc_lo;
    xmax_bc = bc_hi;
}

void
ParticleBoundaries::SetBoundsY (ParticleBoundaryType bc_lo, ParticleBoundaryType bc_hi)
{
    ymin_bc = bc_lo;
    ymax_bc = bc_hi;
}

void
ParticleBoundaries::SetBoundsZ (ParticleBoundaryType bc_lo, ParticleBoundaryType bc_hi)
{
    zmin_bc = bc_lo;
    zmax_bc = bc_hi;
}

bool
ParticleBoundaries::CheckAll (ParticleBoundaryType bc)
{
    return (xmin_bc == bc && xmax_bc == bc
#ifdef WARPX_DIM_3D
         && ymin_bc == bc && ymax_bc == bc
#endif
         && zmin_bc == bc && zmax_bc == bc);
}

void
ParticleBoundaries::BuildReflectionModelParsers ()
{
    /*
    reflection_model_xlo_parser = std::make_unique<amrex::Parser>(makeParser(reflection_model_xlo_str, {"v"}));
    reflection_model_xlo = reflection_model_xlo_parser->compile<1>();
    reflection_model_xhi_parser = std::make_unique<amrex::Parser>(makeParser(reflection_model_xhi_str, {"v"}));
    reflection_model_xhi = reflection_model_xhi_parser->compile<1>();

    reflection_model_ylo_parser = std::make_unique<amrex::Parser>(makeParser(reflection_model_ylo_str, {"v"}));
    reflection_model_ylo = reflection_model_ylo_parser->compile<1>();
    reflection_model_yhi_parser = std::make_unique<amrex::Parser>(makeParser(reflection_model_yhi_str, {"v"}));
    reflection_model_yhi = reflection_model_yhi_parser->compile<1>();
    */
    auto reflection_model_zlo_parser = std::make_unique<amrex::Parser>(makeParser("0", {"v"}));
    reflection_model_zlo = reflection_model_zlo_parser->compile<1>();
    auto reflection_model_zhi_parser = std::make_unique<amrex::Parser>(makeParser("0.5", {"v"}));
    reflection_model_zhi = reflection_model_zhi_parser->compile<1>();
}
