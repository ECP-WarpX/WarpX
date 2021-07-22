#ifndef WARPX_FIELD_ACCESSOR_FUNCTORS_H
#define WARPX_FIELD_ACCESSOR_FUNCTORS_H

#include "WarpX.H"

/**
 * \brief Functor that returns the division of the source m_field Array4 value
          by macroparameter, m_parameter value at the respective (i,j,k,ncomp).
 */
struct FieldAccessorMacroscopic
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    FieldAccessorMacroscopic ( amrex::Array4<amrex::Real const> const a_field,
                               amrex::Array4<amrex::Real const> const a_parameter )
        : m_field(a_field), m_parameter(a_parameter) {}

    /**
     * \brief return field value at (i,j,k,ncomp) scaled by (1/m_parameters(i,j,k,ncomp))
     * \param[in] i      index along x of the Array4, m_field and m_parameter.
     * \param[in] j      index along y of the Array4, m_field and m_parameter.
     * \param[in] k      index along z of the Array4, m_field and m_parameter.
     * \param[in] ncomp  index along fourth component of the Array4, containing field-data
                         to be returned after diving with zero-th component
                         of m_paramter.
     *
     * \return           m_field/m_paramter at (i,j,k,ncomp)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int const i, int const j,
                            int const k, int const ncomp) const noexcept
    {
        return ( m_field(i, j, k, ncomp)  /  m_parameter(i, j, k, 0) ) ;
    }
private:
    /** Array4 of the source field to be scaled and returned by the operator() */
    amrex::Array4<amrex::Real const> const m_field;
    /** Array4 of the macroscopic parameter used to divide m_field in the operator() */
    amrex::Array4<amrex::Real const> const m_parameter;
};


#endif
