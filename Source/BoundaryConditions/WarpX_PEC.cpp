#include "BoundaryConditions/WarpX_PEC.H"
#include "WarpX.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

using namespace amrex::literals;

bool
PEC::isAnyBoundaryPEC() {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC) { return true; }
        if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC) { return true; }
    }
    return false;
}

void
PEC::ApplyPECtoEfield (std::array<amrex::MultiFab*, 3> Efield, const int lev,
                       PatchType patch_type, const bool split_pml_field)
{
    auto& warpx = WarpX::GetInstance();
    amrex::Box domain_box = warpx.Geom(lev).Domain();
    if (patch_type == PatchType::coarse) {
        const amrex::IntVect ref_ratio = ( (lev > 0) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1) );
        domain_box.coarsen(ref_ratio);
    }
    const amrex::IntVect domain_lo = domain_box.smallEnd();
    const amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::GpuArray<int, 3> fbndry_lo;
    amrex::GpuArray<int, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    const amrex::IntVect Ex_nodal = Efield[0]->ixType().toIntVect();
    const amrex::IntVect Ey_nodal = Efield[1]->ixType().toIntVect();
    const amrex::IntVect Ez_nodal = Efield[2]->ixType().toIntVect();
    const amrex::IntVect ng_fieldgather = warpx.get_ng_fieldgather();
    // For each Efield multifab, apply PEC boundary condition to ncomponents
    // If not split E-field, the PEC is applied to the regular Efield used in Maxwell's eq.
    // If split_pml_field is true, then PEC is applied to all the split field components of the tangential field.
    const int nComp_x = Efield[0]->nComp();
    const int nComp_y = Efield[1]->nComp();
    const int nComp_z = Efield[2]->nComp();
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

        // loop over cells and update fields
        amrex::ParallelFor(
            tex, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
                PEC::SetEfieldOnPEC(icomp, domain_lo, domain_hi, iv, n,
                                           Ex, Ex_nodal, fbndry_lo, fbndry_hi);
            },
            tey, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
                PEC::SetEfieldOnPEC(icomp, domain_lo, domain_hi, iv, n,
                                           Ey, Ey_nodal, fbndry_lo, fbndry_hi);
            },
            tez, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
                PEC::SetEfieldOnPEC(icomp, domain_lo, domain_hi, iv, n,
                                           Ez, Ez_nodal, fbndry_lo, fbndry_hi);
            }
        );
    }
}


void
PEC::ApplyPECtoBfield (std::array<amrex::MultiFab*, 3> Bfield, const int lev,
                       PatchType patch_type)
{
    auto& warpx = WarpX::GetInstance();
    amrex::Box domain_box = warpx.Geom(lev).Domain();
    if (patch_type == PatchType::coarse) {
        const amrex::IntVect ref_ratio = ( (lev > 0) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1) );
        domain_box.coarsen(ref_ratio);
    }
    const amrex::IntVect domain_lo = domain_box.smallEnd();
    const amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::GpuArray<int, 3> fbndry_lo;
    amrex::GpuArray<int, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    const amrex::IntVect Bx_nodal = Bfield[0]->ixType().toIntVect();
    const amrex::IntVect By_nodal = Bfield[1]->ixType().toIntVect();
    const amrex::IntVect Bz_nodal = Bfield[2]->ixType().toIntVect();
    const amrex::IntVect ng_fieldgather = warpx.get_ng_fieldgather();
    const int nComp_x = Bfield[0]->nComp();
    const int nComp_y = Bfield[1]->nComp();
    const int nComp_z = Bfield[2]->nComp();

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

        // loop over cells and update fields
        amrex::ParallelFor(
            tbx, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
                PEC::SetBfieldOnPEC(icomp, domain_lo, domain_hi, iv, n,
                                     Bx, Bx_nodal, fbndry_lo, fbndry_hi);
            },
            tby, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
                PEC::SetBfieldOnPEC(icomp, domain_lo, domain_hi, iv, n,
                                     By, By_nodal, fbndry_lo, fbndry_hi);
            },
            tbz, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
                PEC::SetBfieldOnPEC(icomp, domain_lo, domain_hi, iv, n,
                                     Bz, Bz_nodal, fbndry_lo, fbndry_hi);
            }
        );
    }
}


/**
 * \brief Sets the rho field value in cells close to and inside a PEC boundary.
 *        The charge density deposited in the guard cells are either reflected
 *        back into the simulation domain (if a reflecting particle
 *        boundary is used), or the opposite charge density is deposited
 *        back in the domain to capture the effect of an image charge.
 *        The charge density on the PEC boundary is set to 0 while values
 *        in the guard cells are set equal and opposite to their mirror
 *        location inside the domain - representing image charges.
 **/
void
PEC::ApplyPECtoRhofield (amrex::MultiFab* rho, const int lev, PatchType patch_type)
{
    auto& warpx = WarpX::GetInstance();

    amrex::Box domain_box = warpx.Geom(lev).Domain();
    if (patch_type == PatchType::coarse) {
        const amrex::IntVect ref_ratio = ( (lev > 0) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1) );
        domain_box.coarsen(ref_ratio);
    }
    domain_box.convert(rho->ixType());

    amrex::IntVect domain_lo = domain_box.smallEnd();
    amrex::IntVect domain_hi = domain_box.bigEnd();

    amrex::IntVect rho_nodal = rho->ixType().toIntVect();
    amrex::IntVect ng_fieldgather = rho->nGrowVect();

    // Create a copy of the domain which will be extended to include guard
    // cells for boundaries that are NOT PEC
    amrex::Box grown_domain_box = domain_box;

    amrex::GpuArray<GpuArray<bool,2>, AMREX_SPACEDIM> is_pec;
    amrex::GpuArray<bool, AMREX_SPACEDIM> is_tangent_to_bndy;
    amrex::GpuArray<GpuArray<amrex::Real,2>, AMREX_SPACEDIM> psign;
    amrex::GpuArray<GpuArray<int,2>, AMREX_SPACEDIM> mirrorfac;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        is_pec[idim][0] = WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC;
        is_pec[idim][1] = WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC;
        if (!is_pec[idim][0]) { grown_domain_box.growLo(idim, ng_fieldgather[idim]); }
        if (!is_pec[idim][1]) { grown_domain_box.growHi(idim, ng_fieldgather[idim]); }

        // rho values inside guard cells are updated the same as tangential
        // components of the current density
        is_tangent_to_bndy[idim] = true;

        psign[idim][0] = (WarpX::particle_boundary_lo[idim] == ParticleBoundaryType::Reflecting)
                         ? 1._rt : -1._rt;
        psign[idim][1] = (WarpX::particle_boundary_hi[idim] == ParticleBoundaryType::Reflecting)
                         ? 1._rt : -1._rt;
        mirrorfac[idim][0] = 2*domain_lo[idim] - (1 - rho_nodal[idim]);
        mirrorfac[idim][1] = 2*domain_hi[idim] + (1 - rho_nodal[idim]);
    }
    const int nComp = rho->nComp();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*rho); mfi.isValid(); ++mfi) {

        // Get the multifab box including ghost cells
        Box const& fabbox = mfi.fabbox();

        // If grown_domain_box contains fabbox it means there are no PEC
        // boundaries to handle so continue to next box
        if (grown_domain_box.contains(fabbox)) { continue; }

        // Extract field data
        auto const& rho_array = rho->array(mfi);

        // Loop over valid cells (i.e. cells inside the domain)
        amrex::ParallelFor(mfi.validbox(), nComp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            amrex::ignore_unused(k);
#elif (defined WARPX_DIM_1D_Z)
            amrex::ignore_unused(j,k);
#endif
            // Store the array index
            const amrex::IntVect iv(AMREX_D_DECL(i,j,k));

            PEC::SetRhoOrJfieldFromPEC(
                n, iv, rho_array, mirrorfac, psign, is_pec,
                is_tangent_to_bndy, fabbox
            );
        });
    }
}


void
PEC::ApplyPECtoJfield(amrex::MultiFab* Jx, amrex::MultiFab* Jy,
                      amrex::MultiFab* Jz, const int lev,
                      PatchType patch_type)
{
    auto& warpx = WarpX::GetInstance();

    amrex::Box domain_box = warpx.Geom(lev).Domain();
    if (patch_type == PatchType::coarse) {
        const amrex::IntVect ref_ratio = ( (lev > 0) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1) );
        domain_box.coarsen(ref_ratio);
    }

    // Note: force domain box to be nodal to simplify the mirror cell
    // calculations below
    domain_box.convert(IntVect::TheNodeVector());

    amrex::IntVect domain_lo = domain_box.smallEnd();
    amrex::IntVect domain_hi = domain_box.bigEnd();

    // Get the nodal flag for each current component since it is needed to
    // determine the appropriate mirrorfac values below
    amrex::IntVect Jx_nodal = Jx->ixType().toIntVect();
    amrex::IntVect Jy_nodal = Jy->ixType().toIntVect();
    amrex::IntVect Jz_nodal = Jz->ixType().toIntVect();

    // Create a copy of the domain for each component of J which will
    // be extended to include guard cells for boundaries that are NOT PEC
    amrex::Box grown_domain_box = domain_box;

    // The number of ghost cells used is the same for nodal and cell-centered
    // directions of the current density multifab
    const amrex::IntVect ng_fieldgather = Jx->nGrowVect();

    amrex::GpuArray<GpuArray<bool, 2>, AMREX_SPACEDIM> is_pec;
    amrex::GpuArray<GpuArray<bool, AMREX_SPACEDIM>, 3> is_tangent_to_bndy;
    amrex::GpuArray<GpuArray<GpuArray<amrex::Real, 2>, AMREX_SPACEDIM>, 3> psign;
    amrex::GpuArray<GpuArray<GpuArray<int, 2>, AMREX_SPACEDIM>, 3> mirrorfac;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        is_pec[idim][0] = WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC;
        is_pec[idim][1] = WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC;
        if (!is_pec[idim][0]) { grown_domain_box.growLo(idim, ng_fieldgather[idim]); }
        if (!is_pec[idim][1]) { grown_domain_box.growHi(idim, ng_fieldgather[idim]); }

        for (int icomp=0; icomp < 3; ++icomp) {
            // Set the psign value correctly for each current component for each
            // simulation direction
#if (defined WARPX_DIM_1D_Z)
            // For 1D : icomp=0 and icomp=1 (Ex and Ey are tangential to the z boundary)
            //          The logic below ensures that the flags are set right for 1D
            is_tangent_to_bndy[icomp][idim] = (icomp != (idim+2));
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            // For 2D : for icomp==1, (Ey in XZ, Etheta in RZ),
            //          icomp=1 is tangential to both x and z boundaries
            //          The logic below ensures that the flags are set right for 2D
            is_tangent_to_bndy[icomp][idim] = (icomp != AMREX_SPACEDIM*idim);
#else
            is_tangent_to_bndy[icomp][idim] = (icomp != idim);
#endif

            if (is_tangent_to_bndy[icomp][idim]){
                psign[icomp][idim][0] = (WarpX::particle_boundary_lo[idim] == ParticleBoundaryType::Reflecting)
                                        ? 1._rt : -1._rt;
                psign[icomp][idim][1] = (WarpX::particle_boundary_hi[idim] == ParticleBoundaryType::Reflecting)
                                        ? 1._rt : -1._rt;
            }
            else {
                psign[icomp][idim][0] = (WarpX::particle_boundary_lo[idim] == ParticleBoundaryType::Reflecting)
                                        ? -1._rt : 1._rt;
                psign[icomp][idim][1] = (WarpX::particle_boundary_hi[idim] == ParticleBoundaryType::Reflecting)
                                        ? -1._rt : 1._rt;
            }
        }
        // Set the correct mirror cell calculation for each current
        // component. Seeing as domain was set to be nodal, nodal fields have
        // boundaries on domain_lo and domain_hi, but cell-centered fields
        // have valid cells ranging from domain_lo to domain_hi - 1.
        mirrorfac[0][idim][0] = 2*domain_lo[idim] - (1 - Jx_nodal[idim]);
        mirrorfac[0][idim][1] = 2*domain_hi[idim] - (1 - Jx_nodal[idim]);
        mirrorfac[1][idim][0] = 2*domain_lo[idim] - (1 - Jy_nodal[idim]);
        mirrorfac[1][idim][1] = 2*domain_hi[idim] - (1 - Jy_nodal[idim]);
        mirrorfac[2][idim][0] = 2*domain_lo[idim] - (1 - Jz_nodal[idim]);
        mirrorfac[2][idim][1] = 2*domain_hi[idim] - (1 - Jz_nodal[idim]);
    }

    // Each current component is handled separately below, starting with Jx.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Jx); mfi.isValid(); ++mfi) {

        // Get the multifab box including ghost cells
        Box const& fabbox = mfi.fabbox();

        // If grown_domain_box contains fabbox it means there are no PEC
        // boundaries to handle so continue to next box
        grown_domain_box.convert(Jx_nodal);
        if (grown_domain_box.contains(fabbox)) { continue; }

        // Extract field data
        auto const& Jx_array = Jx->array(mfi);

        // Loop over valid cells (i.e. cells inside the domain)
        amrex::ParallelFor(mfi.validbox(), Jx->nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            amrex::ignore_unused(k);
#elif (defined WARPX_DIM_1D_Z)
            amrex::ignore_unused(j,k);
#endif
            // Store the array index
            const amrex::IntVect iv(AMREX_D_DECL(i,j,k));

            PEC::SetRhoOrJfieldFromPEC(
                n, iv, Jx_array, mirrorfac[0], psign[0], is_pec,
                is_tangent_to_bndy[0], fabbox
            );
        });
    }

    // Handle Jy.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Jy); mfi.isValid(); ++mfi) {

        // Get the multifab box including ghost cells
        Box const& fabbox = mfi.fabbox();

        // If grown_domain_box contains fabbox it means there are no PEC
        // boundaries to handle so continue to next box
        grown_domain_box.convert(Jy_nodal);
        if (grown_domain_box.contains(fabbox)) { continue; }

        // Extract field data
        auto const& Jy_array = Jy->array(mfi);

        // Loop over valid cells (i.e. cells inside the domain)
        amrex::ParallelFor(mfi.validbox(), Jy->nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            amrex::ignore_unused(k);
#elif (defined WARPX_DIM_1D_Z)
            amrex::ignore_unused(j,k);
#endif
            // Store the array index
            const amrex::IntVect iv(AMREX_D_DECL(i,j,k));

            PEC::SetRhoOrJfieldFromPEC(
                n, iv, Jy_array, mirrorfac[1], psign[1], is_pec,
                is_tangent_to_bndy[1], fabbox
            );
        });
    }

    // Handle Jz.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Jz); mfi.isValid(); ++mfi) {

        // Get the multifab box including ghost cells
        Box const& fabbox = mfi.fabbox();

        // If grown_domain_box contains fabbox it means there are no PEC
        // boundaries to handle so continue to next box
        grown_domain_box.convert(Jz_nodal);
        if (grown_domain_box.contains(fabbox)) { continue; }

        // Extract field data
        auto const& Jz_array = Jz->array(mfi);

        // Loop over valid cells (i.e. cells inside the domain)
        amrex::ParallelFor(mfi.validbox(), Jz->nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            amrex::ignore_unused(k);
#elif (defined WARPX_DIM_1D_Z)
            amrex::ignore_unused(j,k);
#endif
            // Store the array index
            const amrex::IntVect iv(AMREX_D_DECL(i,j,k));

            PEC::SetRhoOrJfieldFromPEC(
                n, iv, Jz_array, mirrorfac[2], psign[2], is_pec,
                is_tangent_to_bndy[2], fabbox
            );
        });
    }
}

void
PEC::ApplyPECtoElectronPressure (amrex::MultiFab* Pefield, const int lev,
                                 PatchType patch_type)
{
    auto& warpx = WarpX::GetInstance();

    amrex::Box domain_box = warpx.Geom(lev).Domain();
    if (patch_type == PatchType::coarse) {
        const amrex::IntVect ref_ratio = ( (lev > 0) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1) );
        domain_box.coarsen(ref_ratio);
    }
    domain_box.convert(Pefield->ixType());

    amrex::IntVect domain_lo = domain_box.smallEnd();
    amrex::IntVect domain_hi = domain_box.bigEnd();

    amrex::IntVect Pe_nodal = Pefield->ixType().toIntVect();
    amrex::IntVect ng_fieldgather = Pefield->nGrowVect();

    // Create a copy of the domain which will be extended to include guard
    // cells for boundaries that are NOT PEC
    amrex::Box grown_domain_box = domain_box;

    amrex::GpuArray<GpuArray<bool,2>, AMREX_SPACEDIM> is_pec;
    amrex::GpuArray<GpuArray<int,2>, AMREX_SPACEDIM> mirrorfac;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        is_pec[idim][0] = WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC;
        is_pec[idim][1] = WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC;
        if (!is_pec[idim][0]) { grown_domain_box.growLo(idim, ng_fieldgather[idim]); }
        if (!is_pec[idim][1]) { grown_domain_box.growHi(idim, ng_fieldgather[idim]); }

        mirrorfac[idim][0] = 2*domain_lo[idim] - (1 - Pe_nodal[idim]);
        mirrorfac[idim][1] = 2*domain_hi[idim] + (1 - Pe_nodal[idim]);
    }
    const int nComp = Pefield->nComp();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Pefield); mfi.isValid(); ++mfi) {

        // Get the multifab box including ghost cells
        Box const& fabbox = mfi.fabbox();

        // If grown_domain_box contains fabbox it means there are no PEC
        // boundaries to handle so continue to next box
        if (grown_domain_box.contains(fabbox)) { continue; }

        // Extract field data
        auto const& Pe_array = Pefield->array(mfi);

        // Loop over valid cells (i.e. cells inside the domain)
        amrex::ParallelFor(mfi.validbox(), nComp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            amrex::ignore_unused(k);
#elif (defined WARPX_DIM_1D_Z)
            amrex::ignore_unused(j,k);
#endif
            // Store the array index
            const amrex::IntVect iv(AMREX_D_DECL(i,j,k));

            PEC::SetNeumannOnPEC(n, iv, Pe_array, mirrorfac, is_pec, fabbox);
        });
    }
}
