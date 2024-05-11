#ifndef PEC_INSULATOR_H_
#define PEC_INSULATOR_H_

#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

class PEC_Insulator
{
public:

    PEC_Insulator();

    /**
     * \brief The electric field at the insulator boundary is unchanged.
     *        The guard cell values are extrapolated from the valid cell field values.
     *        Outside of the insulator, PEC boundaries are applied.
     *
     * \param[in,out] Efield
     * \param[in]     ng_fieldgather      number of guard cells used by field gather
     * \param[in]     geom                geometry object of level "lev"
     * \param[in]     lev                 level of the Multifab
     * \param[in]     patch_type          coarse or fine
     * \param[in]     ref_ratios          vector containing the refinement ratios of the refinement levels
     * \param[in]     split_pml_field     whether pml the multifab is the regular Efield or
     *                                        split pml field
     */
    void ApplyPEC_InsulatortoEfield (std::array<amrex::MultiFab*, 3> Efield,
                                     amrex::Vector<FieldBoundaryType> const & field_boundary_lo,
                                     amrex::Vector<FieldBoundaryType> const & field_boundary_hi,
                                     amrex::IntVect const & ng_fieldgather, amrex::Geometry const & geom,
                                     int lev, PatchType patch_type, amrex::Vector<amrex::IntVect> const & ref_ratios,
                                     bool split_pml_field = false);
    /**
     * \brief The magnetic field at the insulator boundary are unchanged.
     *        The tangential guard cell values are set to the specified values.
     *        The normal guard cells values are extrapolated from the valid cell field values
     *        Outside of the insulator, PEC boundaries are applied.
     *
     * \param[in,out] Bfield
     * \param[in]     ng_fieldgather      number of guard cells used by field gather
     * \param[in]     geom                geometry object of level "lev"
     * \param[in]     lev                 level of the Multifab
     * \param[in]     patch_type          coarse or fine
     * \param[in]     ref_ratios          vector containing the refinement ratios of the refinement levels
     * \param[in]     time                current time of the simulation
     */
    void ApplyPEC_InsulatortoBfield (std::array<amrex::MultiFab*, 3> Bfield,
                                     amrex::Vector<FieldBoundaryType> const & field_boundary_lo,
                                     amrex::Vector<FieldBoundaryType> const & field_boundary_hi,
                                     amrex::IntVect const & ng_fieldgather, amrex::Geometry const & geom,
                                     int lev, PatchType patch_type, amrex::Vector<amrex::IntVect> const & ref_ratios,
                                     amrex::Real time);

private:

    std::vector<std::unique_ptr<amrex::Parser>> m_insulator_area_lo;
    std::vector<std::unique_ptr<amrex::Parser>> m_insulator_area_hi;

    std::unique_ptr<amrex::Parser> m_By_x_lo, m_Bz_x_lo;
    std::unique_ptr<amrex::Parser> m_By_x_hi, m_Bz_x_hi;
    std::unique_ptr<amrex::Parser> m_Bx_y_lo, m_Bz_y_lo;
    std::unique_ptr<amrex::Parser> m_Bx_y_hi, m_Bz_y_hi;
    std::unique_ptr<amrex::Parser> m_Bx_z_lo, m_By_z_lo;
    std::unique_ptr<amrex::Parser> m_Bx_z_hi, m_By_z_hi;

};
#endif // PEC_INSULATOR_H_
