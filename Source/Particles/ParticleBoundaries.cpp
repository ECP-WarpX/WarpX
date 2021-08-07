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

