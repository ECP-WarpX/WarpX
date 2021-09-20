/* Copyright 2020 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FIELD_ACCESSOR_FUNCTORS_H
#define WARPX_FIELD_ACCESSOR_FUNCTORS_H

#include "WarpX.H"
#include "Utils/CoarsenIO.H"
#include "Utils/WarpXUtil.H"
#include "FieldSolver/FiniteDifferenceSolver/MacroscopicProperties/MacroscopicProperties.H"
#include <AMReX_Array.H>

/**
 * \brief Functor that returns the division of the source m_field Array4 value
          by macroparameter obtained using functor, m_getParameter,
          at the respective (i,j,k,ncomp).
 */
template< typename T_GetMacroparameter>
struct FieldAccessorMacroscopic
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    FieldAccessorMacroscopic ( amrex::Array4<amrex::Real const> const a_field,
                               T_GetMacroparameter const& a_getParameter,
                               amrex::GpuArray<int,3> const& a_field_stag,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const a_domain_lo,
                               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const a_dx )
        : m_field(a_field), m_getParameter(a_getParameter), m_field_stag(a_field_stag),
          m_domain_lo(a_domain_lo), m_dx(a_dx) {}

    /**
     * \brief return field value at (i,j,k,ncomp) scaled by (1/m_getParameter(x,y,z))
     *
     * \param[in] i      index along x of the Array4, m_field and m_parameter.
     * \param[in] j      index along y of the Array4, m_field and m_parameter.
     * \param[in] k      index along z of the Array4, m_field and m_parameter.
     * \param[in] ncomp  index along fourth component of the Array4, containing field-data
                         to be returned after dividing by the macroparameter.
     *
     * \return           m_field/m_getParameter(x,y,z) at (i,j,k,ncomp)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int const i, int const j,
                            int const k, int const ncomp) const noexcept
    {
        amrex::Real x, y, z;
        WarpXUtilAlgo::getCellCoordinates(i, j, k, m_field_stag, m_domain_lo, m_dx, x, y, z);
        return ( m_field(i, j, k, ncomp)  /  m_getParameter(x, y, z) ) ;
    }
private:
    /** Array4 of the source field to be scaled and returned by the operator() */
    amrex::Array4<amrex::Real const> const m_field;
    /** Functor object to return the macroparameter at a given position on the grid.*/
    T_GetMacroparameter const m_getParameter;
    /** Staggering of the field multifab, m_field */
    amrex::GpuArray<int,3> const m_field_stag;
    /** Lower physical coordinates of the simulation domain. */
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const m_domain_lo;
    /** Cell-size array */
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const m_dx;
};


#endif
