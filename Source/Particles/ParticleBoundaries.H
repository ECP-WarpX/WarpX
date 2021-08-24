/* Copyright 2021 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARIES_H_
#define PARTICLEBOUNDARIES_H_

#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_AmrCore.H>

/**
 *  This allows a different boundary condition on each edge of the grid.
 */
struct ParticleBoundaries
{

    ParticleBoundaries () noexcept;

    void Set_reflect_all_velocities (bool flag);
    void SetAll (ParticleBoundaryType bc);

    void SetBoundsX (ParticleBoundaryType bc_lo, ParticleBoundaryType bc_hi);
    void SetBoundsY (ParticleBoundaryType bc_lo, ParticleBoundaryType bc_hi);
    void SetBoundsZ (ParticleBoundaryType bc_lo, ParticleBoundaryType bc_hi);

    bool CheckAll (ParticleBoundaryType bc);

    ParticleBoundaryType xmin_bc;
    ParticleBoundaryType xmax_bc;
    ParticleBoundaryType ymin_bc;
    ParticleBoundaryType ymax_bc;
    ParticleBoundaryType zmin_bc;
    ParticleBoundaryType zmax_bc;

    bool reflect_all_velocities;

};
#endif /*PARTICLEBOUNDARIES_H_*/
