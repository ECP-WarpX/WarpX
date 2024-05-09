#include "BoundaryConditions/PEC_Insulator.H"
#include "Utils/Parser/ParserUtils.H"
#include "WarpX.H"

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

namespace
{
    /**
     * \brief Calculates the number of grid points the given index is beyond the
     *        domain boundary i.e. a value of +1 means the current cell is
     *        outside of the simulation domain by 1 cell. Note that the high
     *        side domain boundary is between cell dom_hi and dom_hi+1 for cell
     *        centered grids and on cell dom_hi+1 for nodal grid. This is why
     *        (dom_hi[idim] + is_nodal[idim]) is used below.
     *
     * \param[in] dom_lo, dom_hi  Domain boundaries
     * \param[in] ijk_vec         Cell coordinates
     * \param[in] is_nodal        Whether the field of interest is nodal
     * \param[in] idim            Dimension of interest
     * \param[in] iside           -1 for low and +1 for high
     *
     * \returns number of grid points beyond the boundary
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int get_cell_count_to_boundary (const amrex::IntVect& dom_lo,
        const amrex::IntVect& dom_hi, const amrex::IntVect& ijk_vec,
        const amrex::IntVect& is_nodal, const int idim, const int iside)
    {
        return ((iside == -1) ? (dom_lo[idim] - ijk_vec[idim])
                             : (ijk_vec[idim] - (dom_hi[idim] + is_nodal[idim])));
    }


    /**
     * \brief The Efield in the guard cells outside the domain boundary are set to
     *        a value extrapolated from the valid cells.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *
     * \param[in] icomp        component of the Efield being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] Efield       field data to be updated if (ijk) is at the boundary or a guard cell
     * \param[in] is_nodal     staggering of the field data being updated.
     * \param[in] is_insulator_lo Specifies whether lower boundaries are insulators
     * \param[in] is_insulator_hi Specifies whether upper boundaries are insulators
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetEfieldOnPEC_Insulator (const int icomp, const amrex::IntVect & dom_lo,
                               const amrex::IntVect &dom_hi,
                               const amrex::IntVect &ijk_vec, const int n,
                               amrex::Array4<amrex::Real> const& Efield,
                               const amrex::IntVect& is_nodal,
                               const amrex::IntVect & is_insulator_lo,
                               const amrex::IntVect & is_insulator_hi,
                               amrex::GpuArray<FieldBoundaryType, 3> fbndry_lo,
                               amrex::GpuArray<FieldBoundaryType, 3> fbndry_hi)
    {
        using namespace amrex::literals;
        amrex::IntVect ijk_mirror = ijk_vec;
        amrex::IntVect ijk_mirrorp1 = ijk_vec;
        bool OnPECBoundary = false;
        bool GuardCell = false;
        bool isInsulatorBoundary = false;
        amrex::Real sign = 1._rt;
        // Loop over all dimensions
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Loop over sides, iside = -1 (lo), iside = +1 (hi)
            for (int iside = -1; iside <= +1; iside += 2) {
                const bool isPEC_InsulatorBoundary = ( (iside == -1)
                    ? fbndry_lo[idim] == FieldBoundaryType::PEC_Insulator
                    : fbndry_hi[idim] == FieldBoundaryType::PEC_Insulator );
                if (isPEC_InsulatorBoundary) {
                    isInsulatorBoundary = ( (iside == -1)
                        ? is_insulator_lo[idim] == 1
                        : is_insulator_hi[idim] == 1 );
                }
                if (isPEC_InsulatorBoundary) {
                    // grid point ijk_vec is ig number of points pass the
                    // domain boundary in direction, idim
                    const int ig = ::get_cell_count_to_boundary(
                        dom_lo, dom_hi, ijk_vec, is_nodal, idim, iside);

#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                    // For 2D : for icomp==1, (Ey in XZ, Etheta in RZ),
                    //          icomp=1 is tangential to both x and z boundaries
                    //          The logic below ensures that the flags are set right for 2D
                    const bool is_tangent_to_boundary = (icomp != 2*idim);
#elif (defined WARPX_DIM_1D_Z)
                    // For 1D : icomp=0 and icomp=1 (Ex and Ey are tangential to the z boundary)
                    //          The logic below ensures that the flags are set right for 1D
                    const bool is_tangent_to_boundary = (icomp != 2);
#else
                    const bool is_tangent_to_boundary = (icomp != idim);
#endif

                    if (ig == 0) {
                        if (is_tangent_to_boundary && is_nodal[idim] == 1) {
                            OnPECBoundary = true;
                        }
                    } else if (ig > 0) {
                        GuardCell = true;
                        // Find mirror location across boundary
                        ijk_mirror[idim] = ( ( iside == -1)
                                        ? (dom_lo[idim] + ig - (1 - is_nodal[idim]))
                                        : (dom_hi[idim] + 1 - ig));
                        ijk_mirrorp1[idim] = 2*ijk_mirror[idim] - ijk_vec[idim];
                        // tangential components are inverted across boundary
                        if (is_tangent_to_boundary) { sign *= -1._rt; }

#if (defined WARPX_DIM_RZ)
                        if (icomp == 0 && idim == 0 && iside == +1) {
                            // Add radial scale so that drEr/dr = 0.
                            // This only works for the first guard cell and with
                            // Er cell centered in r.
                            const amrex::Real rguard = ijk_vec[idim] + 0.5_rt*(1._rt - is_nodal[idim]);
                            const amrex::Real rmirror = ijk_mirror[idim] + 0.5_rt*(1._rt - is_nodal[idim]);
                            // Calculate radial scale factor
                            sign *= rmirror/rguard;
                        }
#endif
                    }
                } // is pec_insulator boundary
            } // loop over iside
        } // loop over dimensions
        if (isInsulatorBoundary) {
            // The value on the boundary is left unmodified
            // The values in the guard cells are extrapolated
            if (GuardCell) {
                Efield(ijk_vec, n) = 2._rt*Efield(ijk_mirror, n) - Efield(ijk_mirrorp1, n);
            }
        } else {
            if (OnPECBoundary) {
                // if ijk_vec is on a PEC boundary in any direction, set Etangential to 0.
                Efield(ijk_vec,n) = 0._rt;
            } else if (GuardCell) {
                Efield(ijk_vec,n) = sign * Efield(ijk_mirror,n);
            }
        }
    }


    /**
     * \brief The nodal Bfield in the guard cells outside the domain boundary are set to
     *        a value extrapolated from the valid cells.
     *        The non-nodel Bfields in the guard cells are set to specified values.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *
     * \param[in] icomp        component of the Bfield being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] Bfield       field data to be updated if (ijk) is at the boundary
                               or a guard cell
     * \param[in] is_nodal     staggering of the field data being updated.
     * \param[in] is_insulator_lo Specifies whether lower boundaries are insulators
     * \param[in] is_insulator_hi Specifies whether upper boundaries are insulators
     * \param[in] B_lo         Specified values of B at the lower boundaries
     * \param[in] B_hi         Specified values of B at the upper boundaries
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetBfieldOnPEC_Insulator (const int icomp, const amrex::IntVect & dom_lo,
                               const amrex::IntVect & dom_hi,
                               const amrex::IntVect & ijk_vec, const int n,
                               amrex::Array4<amrex::Real> const& Bfield,
                               const amrex::IntVect & is_nodal,
                               const amrex::IntVect & is_insulator_lo,
                               const amrex::IntVect & is_insulator_hi,
                               const amrex::RealVect & B_lo,
                               const amrex::RealVect & B_hi,
                               amrex::GpuArray<FieldBoundaryType, 3> fbndry_lo,
                               amrex::GpuArray<FieldBoundaryType, 3> fbndry_hi)
    {
        using namespace amrex::literals;
        amrex::IntVect ijk_mirror = ijk_vec;
        amrex::IntVect ijk_mirrorp1 = ijk_vec;
        bool OnPECBoundary = false;
        bool GuardCell = false;
        bool isInsulatorBoundary = false;
        amrex::Real sign = 1._rt;
        bool is_normal_to_boundary;
        amrex::Real B_value = 0._rt;
        // Loop over all dimensions
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Loop over sides, iside = -1 (lo), iside = +1 (hi)
            for (int iside = -1; iside <= +1; iside += 2) {
                const bool isPEC_InsulatorBoundary = ( (iside == -1)
                    ? fbndry_lo[idim] == FieldBoundaryType::PEC_Insulator
                    : fbndry_hi[idim] == FieldBoundaryType::PEC_Insulator );
                if (isPEC_InsulatorBoundary) {
                    isInsulatorBoundary = ( (iside == -1)
                        ? is_insulator_lo[idim] == 1
                        : is_insulator_hi[idim] == 1 );
                }
                if (isPEC_InsulatorBoundary) {
                    // grid point ijk_vec is ig number of points pass the
                    // domain boundary in direction, idim
                    const int ig = ::get_cell_count_to_boundary(
                        dom_lo, dom_hi, ijk_vec, is_nodal, idim, iside);

#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                    // For 2D : for icomp==1, (By in XZ, Btheta in RZ),
                    //          icomp=1 is not normal to x or z boundary
                    //          The logic below ensures that the flags are set right for 2D
                    is_normal_to_boundary = (icomp == (2*idim));
#elif (defined WARPX_DIM_1D_Z)
                    // For 1D : icomp=0 and icomp=1 (Bx and By are not normal to the z boundary)
                    //          The logic below ensures that the flags are set right for 1D
                    is_normal_to_boundary = (icomp == 2);
#else
                    is_normal_to_boundary = (icomp == idim);
#endif

                    if (ig == 0) {
                        // Only normal component is set to 0
                        if (is_normal_to_boundary && is_nodal[idim]==1) {
                            OnPECBoundary = true;
                        }
                    } else if ( ig > 0) {
                        GuardCell = true;

                        // Mirror location inside the domain by "ig" number of cells
                        // across boundary in direction, idim, and side, iside
                        ijk_mirror[idim] = ( (iside == -1)
                                        ? (dom_lo[idim] + ig - (1 - is_nodal[idim]))
                                        : (dom_hi[idim] + 1 - ig));
                        ijk_mirrorp1[idim] = 2*ijk_mirror[idim] - ijk_vec[idim];
                        if (is_normal_to_boundary) { sign *= -1._rt; }

                        B_value = ( (iside == -1) ? B_lo[idim] : B_hi[idim] );

#if (defined WARPX_DIM_RZ)
                        if (icomp == 0 && idim == 0 && iside == +1) {
                            // Add radial scale so that drBr/dr = 0.
                            const amrex::Real rguard = ijk_vec[idim] + 0.5_rt*(1._rt - is_nodal[idim]);
                            const amrex::Real rmirror = ijk_mirror[idim] + 0.5_rt*(1._rt - is_nodal[idim]);
                            sign *= rmirror/rguard;
                        }
                        if (isInsulatorBoundary && idim == 0 && iside == +1) {
                            // Add radial scale
                            const amrex::Real rguard = ijk_vec[idim] + 0.5_rt;
                            const amrex::Real rhi = dom_hi[idim];
                            B_value *= rhi/rguard;
                        }
#endif
                    }
                } // if pec_insulator Boundary
            } // loop over sides
        } // loop of dimensions

        if (isInsulatorBoundary) {
            if (GuardCell) {
                if (is_normal_to_boundary) {
                    // The values in the guard cells are extrapolated
                    Bfield(ijk_vec, n) = 2._rt*Bfield(ijk_mirror, n) - Bfield(ijk_mirrorp1, n);
                } else {
                    Bfield(ijk_vec, n) = B_value;
                }
            }
        } else {
            if (OnPECBoundary) {
                // if ijk_vec is on a PEC boundary in any direction, set Bnormal to 0.
                Bfield(ijk_vec,n) = 0._rt;
            } else if (GuardCell) {
                // Bnormal and Btangential is set opposite and equal to the value
                // in the mirror location, respectively.
                Bfield(ijk_vec,n) = sign * Bfield(ijk_mirror,n);
            }
        }
    }
}

PEC_Insulator::PEC_Insulator ()
{

    const amrex::ParmParse pp_insulator("insulator");

#if (AMREX_SPACEDIM > 1)
    std::string str_area_x_lo = "0";
    std::string str_area_x_hi = "0";
    utils::parser::Query_parserString( pp_insulator, "area_x_lo(y,z)", str_area_x_lo);
    utils::parser::Query_parserString( pp_insulator, "area_x_hi(y,z)", str_area_x_hi);
    m_insulator_area_lo.push_back(
        std::make_unique<amrex::Parser>(utils::parser::makeParser(str_area_x_lo, {"y", "z"})));
    m_insulator_area_hi.push_back(
        std::make_unique<amrex::Parser>(utils::parser::makeParser(str_area_x_hi, {"y", "z"})));

    std::string str_By_x_lo = "0";
    std::string str_Bz_x_lo = "0";
    std::string str_By_x_hi = "0";
    std::string str_Bz_x_hi = "0";
    utils::parser::Query_parserString(pp_insulator, "By_x_lo(y,z,t)", str_By_x_lo);
    utils::parser::Query_parserString(pp_insulator, "Bz_x_lo(y,z,t)", str_Bz_x_lo);
    utils::parser::Query_parserString(pp_insulator, "By_x_hi(y,z,t)", str_By_x_hi);
    utils::parser::Query_parserString(pp_insulator, "Bz_x_hi(y,z,t)", str_Bz_x_hi);
    m_By_x_lo = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_By_x_lo, {"y", "z", "t"}));
    m_Bz_x_lo = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bz_x_lo, {"y", "z", "t"}));
    m_By_x_hi = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_By_x_hi, {"y", "z", "t"}));
    m_Bz_x_hi = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bz_x_hi, {"y", "z", "t"}));
#endif
#if defined(WARPX_DIM_3D)
    std::string str_area_y_lo = "0";
    std::string str_area_y_hi = "0";
    utils::parser::Query_parserString( pp_insulator, "area_y_lo(x,z)", str_area_y_lo);
    utils::parser::Query_parserString( pp_insulator, "area_y_hi(x,z)", str_area_y_hi);
    m_insulator_area_lo.push_back(
        std::make_unique<amrex::Parser>(utils::parser::makeParser(str_area_y_lo, {"x", "z"})));
    m_insulator_area_hi.push_back(
        std::make_unique<amrex::Parser>(utils::parser::makeParser(str_area_y_hi, {"x", "z"})));

    std::string str_Bx_y_lo = "0";
    std::string str_Bz_y_lo = "0";
    std::string str_Bx_y_hi = "0";
    std::string str_Bz_y_hi = "0";
    utils::parser::Query_parserString(pp_insulator, "Bx_y_lo(x,z,t)", str_Bx_y_lo);
    utils::parser::Query_parserString(pp_insulator, "Bz_y_lo(x,z,t)", str_Bz_y_lo);
    utils::parser::Query_parserString(pp_insulator, "Bx_y_hi(x,z,t)", str_Bx_y_hi);
    utils::parser::Query_parserString(pp_insulator, "Bz_y_hi(x,z,t)", str_Bz_y_hi);
    m_Bx_y_lo = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bx_y_lo, {"x", "z", "t"}));
    m_Bz_y_lo = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bz_y_lo, {"x", "z", "t"}));
    m_Bx_y_hi = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bx_y_hi, {"x", "z", "t"}));
    m_Bz_y_hi = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bz_y_hi, {"x", "z", "t"}));
#endif

    std::string str_area_z_lo = "0";
    std::string str_area_z_hi = "0";
    utils::parser::Query_parserString( pp_insulator, "area_z_lo(x,y)", str_area_z_lo);
    utils::parser::Query_parserString( pp_insulator, "area_z_hi(x,y)", str_area_z_hi);
    m_insulator_area_lo.push_back(
        std::make_unique<amrex::Parser>(utils::parser::makeParser(str_area_z_lo, {"x", "y"})));
    m_insulator_area_hi.push_back(
        std::make_unique<amrex::Parser>(utils::parser::makeParser(str_area_z_hi, {"x", "y"})));

    std::string str_Bx_z_lo = "0";
    std::string str_By_z_lo = "0";
    std::string str_Bx_z_hi = "0";
    std::string str_By_z_hi = "0";
    utils::parser::Query_parserString(pp_insulator, "Bx_z_lo(x,y,t)", str_Bx_z_lo);
    utils::parser::Query_parserString(pp_insulator, "By_z_lo(x,y,t)", str_By_z_lo);
    utils::parser::Query_parserString(pp_insulator, "Bx_z_hi(x,y,t)", str_Bx_z_hi);
    utils::parser::Query_parserString(pp_insulator, "By_z_hi(x,y,t)", str_By_z_hi);
    m_Bx_z_lo = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bx_z_lo, {"x", "y", "t"}));
    m_By_z_lo = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_By_z_lo, {"x", "y", "t"}));
    m_Bx_z_hi = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_Bx_z_hi, {"x", "y", "t"}));
    m_By_z_hi = std::make_unique<amrex::Parser>(utils::parser::makeParser(str_By_z_hi, {"x", "y", "t"}));

}

void
PEC_Insulator::ApplyPEC_InsulatortoEfield (
    std::array<amrex::MultiFab*, 3> Efield,
    const amrex::Vector<FieldBoundaryType>& field_boundary_lo,
    const amrex::Vector<FieldBoundaryType>& field_boundary_hi,
    const amrex::IntVect& ng_fieldgather, const amrex::Geometry& geom,
    const int lev, PatchType patch_type, const amrex::Vector<amrex::IntVect>& ref_ratios,
    const bool split_pml_field)
{
    using namespace amrex::literals;
    amrex::Box domain_box = geom.Domain();
    if (patch_type == PatchType::coarse && (lev > 0)) {
        domain_box.coarsen(ref_ratios[lev-1]);
    }
    const amrex::IntVect domain_lo = domain_box.smallEnd();
    const amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_lo;
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = field_boundary_lo[idim];
        fbndry_hi[idim] = field_boundary_hi[idim];
    }

#if (AMREX_SPACEDIM > 1)
    amrex::ParserExecutor<2> area_parsers_x_lo = m_insulator_area_lo[0]->compile<2>();
    amrex::ParserExecutor<2> area_parsers_x_hi = m_insulator_area_hi[0]->compile<2>();
#endif
#if defined(WARPX_DIM_3D)
    amrex::ParserExecutor<2> area_parsers_y_lo = m_insulator_area_lo[1]->compile<2>();
    amrex::ParserExecutor<2> area_parsers_y_hi = m_insulator_area_hi[1]->compile<2>();
#endif
    amrex::ParserExecutor<2> area_parsers_z_lo = m_insulator_area_lo[WARPX_ZINDEX]->compile<2>();
    amrex::ParserExecutor<2> area_parsers_z_hi = m_insulator_area_hi[WARPX_ZINDEX]->compile<2>();

    const amrex::IntVect Ex_nodal = Efield[0]->ixType().toIntVect();
    const amrex::IntVect Ey_nodal = Efield[1]->ixType().toIntVect();
    const amrex::IntVect Ez_nodal = Efield[2]->ixType().toIntVect();
    // For each Efield multifab, apply boundary condition to ncomponents
    // If not split E-field, the boundary condition is applied to the regular Efield used in Maxwell's eq.
    // If split_pml_field is true, then boundary condition is applied to all the split field components of the tangential field.
    const int nComp_x = Efield[0]->nComp();
    const int nComp_y = Efield[1]->nComp();
    const int nComp_z = Efield[2]->nComp();

    const std::array<amrex::Real,3>& dx = WarpX::CellSize(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Efield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Extract field data
        amrex::Array4<amrex::Real> const& Ex = Efield[0]->array(mfi);
        amrex::Array4<amrex::Real> const& Ey = Efield[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Ez = Efield[2]->array(mfi);

        // Extract tileboxes for which to loop
        // if split field, the box includes nodal flag
        // For E-field used in Maxwell's update, nodal flag plus cells that particles
        // gather fields from in the guard-cell region are included.
        // Note that for simulations without particles or laser, ng_field_gather is 0
        // and the guard-cell values of the E-field multifab will not be modified.
        amrex::Box const& tex = (split_pml_field) ? mfi.tilebox(Efield[0]->ixType().toIntVect())
                                                  : mfi.tilebox(Efield[0]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tey = (split_pml_field) ? mfi.tilebox(Efield[1]->ixType().toIntVect())
                                                  : mfi.tilebox(Efield[1]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tez = (split_pml_field) ? mfi.tilebox(Efield[2]->ixType().toIntVect())
                                                  : mfi.tilebox(Efield[2]->ixType().toIntVect(), ng_fieldgather);

        std::array<amrex::Real, 3> const& xyzmin_x = WarpX::LowerCorner(tex, lev, 0._rt);
        std::array<amrex::Real, 3> const& xyzmin_y = WarpX::LowerCorner(tey, lev, 0._rt);
        std::array<amrex::Real, 3> const& xyzmin_z = WarpX::LowerCorner(tez, lev, 0._rt);
        amrex::IntVect const lo_x = tex.smallEnd();
        amrex::IntVect const lo_y = tey.smallEnd();
        amrex::IntVect const lo_z = tez.smallEnd();

        // loop over cells and update fields
        amrex::ParallelFor(
            tex, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j, k);
#endif

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real const x = (AMREX_SPACEDIM > 1 ? xyzmin_x[0] + (iv[0] - lo_x[0])*dx[0] : 0.);
                amrex::Real const y = (AMREX_SPACEDIM == 3 ? xyzmin_x[1] + (iv[1] - lo_x[1])*dx[1] : 0.);
#if (AMREX_SPACEDIM > 1)
                amrex::Real const z = xyzmin_x[2] + (iv[WARPX_ZINDEX] - lo_x[WARPX_ZINDEX])*dx[2];
#endif

                amrex::IntVect is_insulator_lo;
                amrex::IntVect is_insulator_hi;
#if (AMREX_SPACEDIM > 1)
                is_insulator_lo[0] = (area_parsers_x_lo(y, z) > 0.);
                is_insulator_hi[0] = (area_parsers_x_hi(y, z) > 0.);
#endif
#if defined(WARPX_DIM_3D)
                is_insulator_lo[1] = (area_parsers_y_lo(x, z) > 0.);
                is_insulator_hi[1] = (area_parsers_y_hi(x, z) > 0.);
#endif
                is_insulator_lo[WARPX_ZINDEX] = (area_parsers_z_lo(x, y) > 0.);
                is_insulator_hi[WARPX_ZINDEX] = (area_parsers_z_hi(x, y) > 0.);

                const int icomp = 0;
                ::SetEfieldOnPEC_Insulator(icomp, domain_lo, domain_hi, iv, n,
                                       Ex, Ex_nodal, is_insulator_lo, is_insulator_hi,
                                       fbndry_lo, fbndry_hi);
            },
            tey, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j, k);
#endif

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real const x = (AMREX_SPACEDIM > 1 ? xyzmin_y[0] + (iv[0] - lo_y[0])*dx[0] : 0.);
                amrex::Real const y = (AMREX_SPACEDIM == 3 ? xyzmin_y[1] + (iv[1] - lo_y[1])*dx[1] : 0.);
#if (AMREX_SPACEDIM > 1)
                amrex::Real const z = xyzmin_y[2] + (iv[WARPX_ZINDEX] - lo_y[WARPX_ZINDEX])*dx[2];
#endif

                amrex::IntVect is_insulator_lo;
                amrex::IntVect is_insulator_hi;
#if (AMREX_SPACEDIM > 1)
                is_insulator_lo[0] = (area_parsers_x_lo(y, z) > 0.);
                is_insulator_hi[0] = (area_parsers_x_hi(y, z) > 0.);
#endif
#if defined(WARPX_DIM_3D)
                is_insulator_lo[1] = (area_parsers_y_lo(x, z) > 0.);
                is_insulator_hi[1] = (area_parsers_y_hi(x, z) > 0.);
#endif
                is_insulator_lo[WARPX_ZINDEX] = (area_parsers_z_lo(x, y) > 0.);
                is_insulator_hi[WARPX_ZINDEX] = (area_parsers_z_hi(x, y) > 0.);

                const int icomp = 1;
                ::SetEfieldOnPEC_Insulator(icomp, domain_lo, domain_hi, iv, n,
                                       Ey, Ey_nodal, is_insulator_lo, is_insulator_hi,
                                       fbndry_lo, fbndry_hi);
            },
            tez, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j, k);
#endif

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real const x = (AMREX_SPACEDIM > 1 ? xyzmin_z[0] + (iv[0] - lo_z[0])*dx[0] : 0.);
                amrex::Real const y = (AMREX_SPACEDIM == 3 ? xyzmin_z[1] + (iv[1] - lo_z[1])*dx[1] : 0.);
#if (AMREX_SPACEDIM > 1)
                amrex::Real const z = xyzmin_z[2] + (iv[WARPX_ZINDEX] - lo_z[WARPX_ZINDEX])*dx[2];
#endif

                amrex::IntVect is_insulator_lo;
                amrex::IntVect is_insulator_hi;
#if (AMREX_SPACEDIM > 1)
                is_insulator_lo[0] = (area_parsers_x_lo(y, z) > 0.);
                is_insulator_hi[0] = (area_parsers_x_hi(y, z) > 0.);
#endif
#if defined(WARPX_DIM_3D)
                is_insulator_lo[1] = (area_parsers_y_lo(x, z) > 0.);
                is_insulator_hi[1] = (area_parsers_y_hi(x, z) > 0.);
#endif
                is_insulator_lo[WARPX_ZINDEX] = (area_parsers_z_lo(x, y) > 0.);
                is_insulator_hi[WARPX_ZINDEX] = (area_parsers_z_hi(x, y) > 0.);

                const int icomp = 2;
                ::SetEfieldOnPEC_Insulator(icomp, domain_lo, domain_hi, iv, n,
                                       Ez, Ez_nodal, is_insulator_lo, is_insulator_hi,
                                       fbndry_lo, fbndry_hi);
            }
        );
    }
}


void
PEC_Insulator::ApplyPEC_InsulatortoBfield (
    std::array<amrex::MultiFab*, 3> Bfield,
    const amrex::Vector<FieldBoundaryType>& field_boundary_lo,
    const amrex::Vector<FieldBoundaryType>& field_boundary_hi,
    const amrex::IntVect& ng_fieldgather, const amrex::Geometry& geom,
    const int lev, PatchType patch_type, const amrex::Vector<amrex::IntVect>& ref_ratios,
    const amrex::Real time)
{
    using namespace amrex::literals;
    amrex::Box domain_box = geom.Domain();
    if (patch_type == PatchType::coarse && (lev > 0)) {
        domain_box.coarsen(ref_ratios[lev-1]);
    }
    const amrex::IntVect domain_lo = domain_box.smallEnd();
    const amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_lo;
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = field_boundary_lo[idim];
        fbndry_hi[idim] = field_boundary_hi[idim];
    }

#if (AMREX_SPACEDIM > 1)
    amrex::ParserExecutor<2> area_parsers_x_lo = m_insulator_area_lo[0]->compile<2>();
    amrex::ParserExecutor<2> area_parsers_x_hi = m_insulator_area_hi[0]->compile<2>();
#endif
#if defined(WARPX_DIM_3D)
    amrex::ParserExecutor<2> area_parsers_y_lo = m_insulator_area_lo[1]->compile<2>();
    amrex::ParserExecutor<2> area_parsers_y_hi = m_insulator_area_hi[1]->compile<2>();
#endif
    amrex::ParserExecutor<2> area_parsers_z_lo = m_insulator_area_lo[WARPX_ZINDEX]->compile<2>();
    amrex::ParserExecutor<2> area_parsers_z_hi = m_insulator_area_hi[WARPX_ZINDEX]->compile<2>();

#if (AMREX_SPACEDIM > 1)
    amrex::ParserExecutor<3> By_x_lo_parser = m_By_x_lo->compile<3>();
    amrex::ParserExecutor<3> Bz_x_lo_parser = m_Bz_x_lo->compile<3>();
    amrex::ParserExecutor<3> By_x_hi_parser = m_By_x_hi->compile<3>();
    amrex::ParserExecutor<3> Bz_x_hi_parser = m_Bz_x_hi->compile<3>();
#endif
#if defined(WARPX_DIM_3D)
    amrex::ParserExecutor<3> Bx_y_lo_parser = m_Bx_y_lo->compile<3>();
    amrex::ParserExecutor<3> Bz_y_lo_parser = m_Bz_y_lo->compile<3>();
    amrex::ParserExecutor<3> Bx_y_hi_parser = m_Bx_y_hi->compile<3>();
    amrex::ParserExecutor<3> Bz_y_hi_parser = m_Bz_y_hi->compile<3>();
#endif
    amrex::ParserExecutor<3> Bx_z_lo_parser = m_Bx_z_lo->compile<3>();
    amrex::ParserExecutor<3> By_z_lo_parser = m_By_z_lo->compile<3>();
    amrex::ParserExecutor<3> Bx_z_hi_parser = m_Bx_z_hi->compile<3>();
    amrex::ParserExecutor<3> By_z_hi_parser = m_By_z_hi->compile<3>();

    const amrex::IntVect Bx_nodal = Bfield[0]->ixType().toIntVect();
    const amrex::IntVect By_nodal = Bfield[1]->ixType().toIntVect();
    const amrex::IntVect Bz_nodal = Bfield[2]->ixType().toIntVect();
    const int nComp_x = Bfield[0]->nComp();
    const int nComp_y = Bfield[1]->nComp();
    const int nComp_z = Bfield[2]->nComp();

    const std::array<amrex::Real,3>& dx = WarpX::CellSize(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Bfield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Extract field data
        amrex::Array4<amrex::Real> const& Bx = Bfield[0]->array(mfi);
        amrex::Array4<amrex::Real> const& By = Bfield[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Bz = Bfield[2]->array(mfi);

        // Extract tileboxes for which to loop
        // For B-field used in Maxwell's update, nodal flag plus cells that particles
        // gather fields from in the guard-cell region are included.
        // Note that for simulations without particles or laser, ng_field_gather is 0
        // and the guard-cell values of the B-field multifab will not be modified.
        amrex::Box const& tbx = mfi.tilebox(Bfield[0]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tby = mfi.tilebox(Bfield[1]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tbz = mfi.tilebox(Bfield[2]->ixType().toIntVect(), ng_fieldgather);

        std::array<amrex::Real, 3> const& xyzmin_x = WarpX::LowerCorner(tbx, lev, 0._rt);
        std::array<amrex::Real, 3> const& xyzmin_y = WarpX::LowerCorner(tby, lev, 0._rt);
        std::array<amrex::Real, 3> const& xyzmin_z = WarpX::LowerCorner(tbz, lev, 0._rt);
        amrex::IntVect const lo_x = tbx.smallEnd();
        amrex::IntVect const lo_y = tby.smallEnd();
        amrex::IntVect const lo_z = tbz.smallEnd();

        // loop over cells and update fields
        amrex::ParallelFor(
            tbx, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j, k);
#endif

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real const x = (AMREX_SPACEDIM > 1 ? xyzmin_x[0] + (iv[0] - lo_x[0])*dx[0] : 0.);
                amrex::Real const y = (AMREX_SPACEDIM == 3 ? xyzmin_x[1] + (iv[1] - lo_x[1])*dx[1] : 0.);
#if (AMREX_SPACEDIM > 1)
                amrex::Real const z = xyzmin_x[2] + (iv[WARPX_ZINDEX] - lo_x[WARPX_ZINDEX])*dx[2];
#endif

                amrex::IntVect is_insulator_lo;
                amrex::IntVect is_insulator_hi;
                amrex::RealVect B_lo, B_hi;
#if (AMREX_SPACEDIM > 1)
                is_insulator_lo[0] = (area_parsers_x_lo(y, z) > 0.);
                is_insulator_hi[0] = (area_parsers_x_hi(y, z) > 0.);
                B_lo[0] = 0.;  // Will be unused
                B_hi[0] = 0.;  // Will be unused
#endif
#if defined(WARPX_DIM_3D)
                is_insulator_lo[1] = (area_parsers_y_lo(x, z) > 0.);
                is_insulator_hi[1] = (area_parsers_y_hi(x, z) > 0.);
                B_lo[1] = Bx_y_lo_parser(x, z, time);
                B_hi[1] = Bx_y_hi_parser(x, z, time);
#endif
                is_insulator_lo[WARPX_ZINDEX] = (area_parsers_z_lo(x, y) > 0.);
                is_insulator_hi[WARPX_ZINDEX] = (area_parsers_z_hi(x, y) > 0.);
                B_lo[WARPX_ZINDEX] = Bx_z_lo_parser(x, y, time);
                B_hi[WARPX_ZINDEX] = Bx_z_hi_parser(x, y, time);

                const int icomp = 0;
                ::SetBfieldOnPEC_Insulator(icomp, domain_lo, domain_hi, iv, n,
                                       Bx, Bx_nodal, is_insulator_lo, is_insulator_hi,
                                       B_lo, B_hi,
                                       fbndry_lo, fbndry_hi);
            },
            tby, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j, k);
#endif

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real const x = (AMREX_SPACEDIM > 1 ? xyzmin_y[0] + (iv[0] - lo_y[0])*dx[0] : 0.);
                amrex::Real const y = (AMREX_SPACEDIM == 3 ? xyzmin_y[1] + (iv[1] - lo_y[1])*dx[1] : 0.);
#if (AMREX_SPACEDIM > 1)
                amrex::Real const z = xyzmin_y[2] + (iv[WARPX_ZINDEX] - lo_y[WARPX_ZINDEX])*dx[2];
#endif

                amrex::IntVect is_insulator_lo;
                amrex::IntVect is_insulator_hi;
                amrex::RealVect B_lo, B_hi;
#if (AMREX_SPACEDIM > 1)
                is_insulator_lo[0] = (area_parsers_x_lo(y, z) > 0.);
                is_insulator_hi[0] = (area_parsers_x_hi(y, z) > 0.);
                B_lo[0] = By_x_lo_parser(y, z, time);
                B_hi[0] = By_x_hi_parser(y, z, time);
#endif
#if defined(WARPX_DIM_3D)
                is_insulator_lo[1] = (area_parsers_y_lo(x, z) > 0.);
                is_insulator_hi[1] = (area_parsers_y_hi(x, z) > 0.);
                B_lo[1] = 0.;  // Will be unused
                B_hi[1] = 0.;  // Will be unused
#endif
                is_insulator_lo[WARPX_ZINDEX] = (area_parsers_z_lo(x, y) > 0.);
                is_insulator_hi[WARPX_ZINDEX] = (area_parsers_z_hi(x, y) > 0.);
                B_lo[WARPX_ZINDEX] = By_z_lo_parser(x, y, time);
                B_hi[WARPX_ZINDEX] = By_z_hi_parser(x, y, time);

                const int icomp = 1;
                ::SetBfieldOnPEC_Insulator(icomp, domain_lo, domain_hi, iv, n,
                                       By, By_nodal, is_insulator_lo, is_insulator_hi,
                                       B_lo, B_hi,
                                       fbndry_lo, fbndry_hi);
            },
            tbz, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j, k);
#endif

                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real const x = (AMREX_SPACEDIM > 1 ? xyzmin_z[0] + (iv[0] - lo_z[0])*dx[0] : 0.);
                amrex::Real const y = (AMREX_SPACEDIM == 3 ? xyzmin_z[1] + (iv[1] - lo_z[1])*dx[1] : 0.);
#if (AMREX_SPACEDIM > 1)
                amrex::Real const z = xyzmin_z[2] + (iv[WARPX_ZINDEX] - lo_z[WARPX_ZINDEX])*dx[2];
#endif

                amrex::IntVect is_insulator_lo;
                amrex::IntVect is_insulator_hi;
                amrex::RealVect B_lo, B_hi;
#if (AMREX_SPACEDIM > 1)
                is_insulator_lo[0] = (area_parsers_x_lo(y, z) > 0.);
                is_insulator_hi[0] = (area_parsers_x_hi(y, z) > 0.);
                B_lo[0] = Bz_x_lo_parser(y, z, time);
                B_hi[0] = Bz_x_hi_parser(y, z, time);
#endif
#if defined(WARPX_DIM_3D)
                is_insulator_lo[1] = (area_parsers_y_lo(x, z) > 0.);
                is_insulator_hi[1] = (area_parsers_y_hi(x, z) > 0.);
                B_lo[1] = Bz_y_lo_parser(x, z, time);
                B_hi[1] = Bz_y_hi_parser(x, z, time);
#endif
                is_insulator_lo[WARPX_ZINDEX] = (area_parsers_z_lo(x, y) > 0.);
                is_insulator_hi[WARPX_ZINDEX] = (area_parsers_z_hi(x, y) > 0.);
                B_lo[WARPX_ZINDEX] = 0.;  // Will be unused
                B_hi[WARPX_ZINDEX] = 0.;  // Will be unused

                const int icomp = 2;
                ::SetBfieldOnPEC_Insulator(icomp, domain_lo, domain_hi, iv, n,
                                       Bz, Bz_nodal, is_insulator_lo, is_insulator_hi,
                                       B_lo, B_hi,
                                       fbndry_lo, fbndry_hi);
            }
        );
    }
}
