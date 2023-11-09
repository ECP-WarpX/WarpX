/* Copyright 2020 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FIELD_ACCESSOR_FUNCTORS_H
#define WARPX_FIELD_ACCESSOR_FUNCTORS_H

#include "WarpX.H"
#include "FieldSolver/FiniteDifferenceSolver/MacroscopicProperties/MacroscopicProperties.H"

#include <AMReX_Array.H>


/**
 * \brief Functor that returns the division of the source m_field Array4 value
          by macroparameter obtained using m_parameter, at the respective (i,j,k).
 */
struct FieldAccessorMacroscopic
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    FieldAccessorMacroscopic ( amrex::Array4<amrex::Real const> const a_field,
                               amrex::Array4<amrex::Real> const& a_parameter)
        : m_field(a_field), m_parameter(a_parameter) {}

    /**
     * \brief return field value at (i,j,k,ncomp) scaled by (1/m_getParameter(x,y,z))
     *
     * \param[in] i      index along x of the Array4, m_field and m_parameter.
     * \param[in] j      index along y of the Array4, m_field and m_parameter.
     * \param[in] k      index along z of the Array4, m_field and m_parameter.
     * \param[in] ncomp  index along fourth component of the Array4, containing field-data
     *                   to be returned after dividing by the macroparameter.
     *
     * \return           m_field/m_getParameter(x,y,z) at (i,j,k,ncomp)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k, int ncomp) const noexcept
    {
        // return ( m_field(i, j, k, ncomp)  /  m_parameter(i, j, k) ) ;
        return (   /  m_parameter(i, j, k) ) * m_field(i, j, k, ncomp);
    }
private:
    /** Array4 of the source field to be scaled and returned by the operator() */
    amrex::Array4<amrex::Real const> const m_field;
    /** Array4 of the macroscopic parameter used to divide m_field in the operator() */
    amrex::Array4<amrex::Real const> const m_parameter;

    int dir;
};

struct HFieldFromBAndM
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k, int ncomp) const noexcept
    {
        return PhysConst::mu0inv * m_B_fc(i,j,k,n) - m_M_cc(i,j,k,n);
    }

private:

    amrex::Array4<const amrex::Real>& m_B_fc;
    amrex::Array4<const amrex::Real>& m_M_cc;


};


#endif
