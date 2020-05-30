/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote, Maxence Thevenet,
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang, Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldIO.H"
#include "WarpX.H"
#include "Utils/CoarsenIO.H"
#include "Utils/WarpXUtil.H"

#ifdef WARPX_USE_PSATD
#   include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif

#include <AMReX_FillPatchUtil_F.H>
#include <AMReX_Interpolater.H>

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif

using namespace amrex;

#ifdef WARPX_USE_OPENPMD
/** \brief For a given field that is to be written to an openPMD file,
 * set the metadata that indicates the physical unit.
 */
void
setOpenPMDUnit( openPMD::Mesh mesh, const std::string field_name )
{
    if (field_name[0] == 'E'){  // Electric field
        mesh.setUnitDimension({
            {openPMD::UnitDimension::L,  1},
            {openPMD::UnitDimension::M,  1},
            {openPMD::UnitDimension::T, -3},
            {openPMD::UnitDimension::I, -1},
        });
    } else if (field_name[0] == 'B'){ // Magnetic field
        mesh.setUnitDimension({
            {openPMD::UnitDimension::M,  1},
            {openPMD::UnitDimension::I, -1},
            {openPMD::UnitDimension::T, -2}
        });
    } else if (field_name[0] == 'j'){ // current
        mesh.setUnitDimension({
            {openPMD::UnitDimension::L, -2},
            {openPMD::UnitDimension::I,  1},
        });
    } else if (field_name.substr(0,3) == "rho"){ // charge density
        mesh.setUnitDimension({
            {openPMD::UnitDimension::L, -3},
            {openPMD::UnitDimension::I,  1},
            {openPMD::UnitDimension::T,  1},
        });
    }
}


/** \brief
 * Convert an IntVect to a std::vector<std::uint64_t>
 * and reverse the order of the elements
 * (used for compatibility with the openPMD API)
 */
std::vector<std::uint64_t>
getReversedVec( const IntVect& v )
{
  // Convert the IntVect v to and std::vector u
  std::vector<std::uint64_t> u = {
    AMREX_D_DECL(
                 static_cast<std::uint64_t>(v[0]),
                 static_cast<std::uint64_t>(v[1]),
                 static_cast<std::uint64_t>(v[2])
                 )
  };
  // Reverse the order of elements, if v corresponds to the indices of a
  // Fortran-order array (like an AMReX FArrayBox)
  // but u is intended to be used with a C-order API (like openPMD)
  std::reverse( u.begin(), u.end() );

  return u;
}

/** \brief
 * Convert Real* pointer to a std::vector<double>,
 * and reverse the order of the elements
 * (used for compatibility with the openPMD API)
 */
std::vector<double>
getReversedVec( const Real* v )
{
  // Convert Real* v to and std::vector u
  std::vector<double> u = {
    AMREX_D_DECL(
                 static_cast<double>(v[0]),
                 static_cast<double>(v[1]),
                 static_cast<double>(v[2])
                 )
  };
  // Reverse the order of elements, if v corresponds to the indices of a
  // Fortran-order array (like an AMReX FArrayBox)
  // but u is intended to be used with a C-order API (like openPMD)
  std::reverse( u.begin(), u.end() );

  return u;
}

#endif // WARPX_USE_OPENPMD

#ifdef WARPX_DIM_RZ
void
ConstructTotalRZVectorField (const std::array< std::unique_ptr<MultiFab>, 3 >& vector_total,
                             const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field)
{
    // Sum over the real components, giving quantity at theta=0
    MultiFab::Copy(*vector_total[0], *vector_field[0], 0, 0, 1, vector_field[0]->nGrowVect());
    MultiFab::Copy(*vector_total[1], *vector_field[1], 0, 0, 1, vector_field[1]->nGrowVect());
    MultiFab::Copy(*vector_total[2], *vector_field[2], 0, 0, 1, vector_field[2]->nGrowVect());
    for (int ic=1 ; ic < vector_field[0]->nComp() ; ic += 2) {
        MultiFab::Add(*vector_total[0], *vector_field[0], ic, 0, 1, vector_field[0]->nGrowVect());
        MultiFab::Add(*vector_total[1], *vector_field[1], ic, 0, 1, vector_field[1]->nGrowVect());
        MultiFab::Add(*vector_total[2], *vector_field[2], ic, 0, 1, vector_field[2]->nGrowVect());
    }
}

void
ConstructTotalRZScalarField (MultiFab& scalar_total,
                            const MultiFab& scalar_field)
{
    // Sum over the real components, giving quantity at theta=0
    MultiFab::Copy(scalar_total, scalar_field, 0, 0, 1, scalar_field.nGrowVect());
    for (int ic=1 ; ic < scalar_field.nComp() ; ic += 2) {
        MultiFab::Add(scalar_total, scalar_field, ic, 0, 1, scalar_field.nGrowVect());
    }
}
#endif

/** \brief Takes an array of 3 MultiFab `vector_field`
 * (representing the x, y, z components of a vector),
 * averages it to the cell center, and stores the
 * resulting MultiFab in mf_avg (in the components dcomp to dcomp+2)
 * Should only be used for BTD now.
 */
void
AverageAndPackVectorField( MultiFab& mf_avg,
                           const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field,
                           const DistributionMapping& dm,
                           const int dcomp, const int ngrow )
{
#ifndef WARPX_DIM_RZ
    (void)dm;
#endif

#ifdef WARPX_DIM_RZ
    // Note that vector_total is declared in the same way as
    // vector_field so that it can be handled the same way.
    std::array<std::unique_ptr<MultiFab>,3> vector_total;
    if (vector_field[0]->nComp() > 1) {
        // With the RZ solver, if there are more than one component, the total
        // fields needs to be constructed in temporary MultiFabs.
        vector_total[0].reset(new MultiFab(vector_field[0]->boxArray(), dm, 1, vector_field[0]->nGrowVect()));
        vector_total[1].reset(new MultiFab(vector_field[1]->boxArray(), dm, 1, vector_field[1]->nGrowVect()));
        vector_total[2].reset(new MultiFab(vector_field[2]->boxArray(), dm, 1, vector_field[2]->nGrowVect()));
        ConstructTotalRZVectorField(vector_total, vector_field);
    } else {
        // Create aliases of the MultiFabs
        vector_total[0].reset(new MultiFab(*vector_field[0], amrex::make_alias, 0, 1));
        vector_total[1].reset(new MultiFab(*vector_field[1], amrex::make_alias, 0, 1));
        vector_total[2].reset(new MultiFab(*vector_field[2], amrex::make_alias, 0, 1));
    }
#else
    const std::array<std::unique_ptr<MultiFab>,3> &vector_total = vector_field;
#endif

    CoarsenIO::Coarsen( mf_avg, *(vector_total[0]), dcomp  , 0, 1, ngrow );
    CoarsenIO::Coarsen( mf_avg, *(vector_total[1]), dcomp+1, 0, 1, ngrow );
    CoarsenIO::Coarsen( mf_avg, *(vector_total[2]), dcomp+2, 0, 1, ngrow );
}

/** \brief Take a MultiFab `scalar_field`
 * averages it to the cell center, and stores the
 * resulting MultiFab in mf_avg (in the components dcomp)
 */
void
AverageAndPackScalarField (MultiFab& mf_avg,
                           const MultiFab & scalar_field,
                           const DistributionMapping& dm,
                           const int dcomp, const int ngrow )
{

#ifdef WARPX_DIM_RZ
    MultiFab *scalar_total;
    if (scalar_field.nComp() > 1) {
        // With the RZ solver, there are more than one component, so the total
        // fields needs to be constructed in temporary a MultiFab.
        scalar_total = new MultiFab(scalar_field.boxArray(), dm, 1, scalar_field.nGrowVect());
        ConstructTotalRZScalarField(*scalar_total, scalar_field);
    } else {
        scalar_total = new MultiFab(scalar_field, amrex::make_alias, 0, 1);
    }
#else
    const MultiFab *scalar_total = &scalar_field;
#endif

    // Check the type of staggering of the 3-component `vector_field`
    // and average accordingly:
    // - Fully cell-centered field (no average needed; simply copy)
    if ( scalar_total->is_cell_centered() ){
        MultiFab::Copy( mf_avg, *scalar_total, 0, dcomp, 1, ngrow);
    } else if ( scalar_total->is_nodal() ){
        // - Fully nodal
        CoarsenIO::Coarsen( mf_avg, *scalar_total, dcomp, 0, 1, ngrow );
    } else {
        amrex::Abort("Unknown staggering.");
    }
}

/** \brief Write the data from MultiFab `F` into the file `filename`
 *  as a raw field (i.e. no interpolation to cell centers).
 *  Write guard cells if `plot_guards` is True.
 */
void
WriteRawField( const MultiFab& F, const DistributionMapping& dm,
                const std::string& filename,
                const std::string& level_prefix,
                const std::string& field_name,
                const int lev, const bool plot_guards )
{
    std::string prefix = amrex::MultiFabFileFullPrefix(lev,
                            filename, level_prefix, field_name);

    if (plot_guards) {
        // Dump original MultiFab F
        VisMF::Write(F, prefix);
    } else {
        // Copy original MultiFab into one that does not have guard cells
        MultiFab tmpF( F.boxArray(), dm, F.nComp(), 0);
        MultiFab::Copy(tmpF, F, 0, 0, F.nComp(), 0);
        VisMF::Write(tmpF, prefix);
    }

}

/** \brief Write a multifab of the same shape as `F` but filled with 0.
 *  (The shape includes guard cells if `plot_guards` is True.)
 *  This is mainly needed because the yt reader requires all levels of the
 *  coarse/fine patch to be written, but WarpX does not have data for
 *  the coarse patch of level 0 (meaningless).
 */
void
WriteZeroRawField( const MultiFab& F, const DistributionMapping& dm,
                const std::string& filename,
                const std::string& level_prefix,
                const std::string& field_name,
                const int lev, const int ng )
{
    std::string prefix = amrex::MultiFabFileFullPrefix(lev,
                            filename, level_prefix, field_name);

    MultiFab tmpF(F.boxArray(), dm, F.nComp(), ng);
    tmpF.setVal(0.);
    VisMF::Write(tmpF, prefix);
}

/** \brief Write the coarse scalar multifab `F_cp` to the file `filename`
 *  *after* sampling/interpolating its value on the fine grid corresponding
 *  to `F_fp`. This is mainly needed because the yt reader requires the
 *  coarse and fine patch to have the same shape.
 */
void
WriteCoarseScalar( const std::string field_name,
    const std::unique_ptr<MultiFab>& F_cp,
    const std::unique_ptr<MultiFab>& F_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards,
    const int icomp )
{
    int ng = 0;
    if (plot_guards) ng = F_fp->nGrow();

    if (lev == 0) {
        // No coarse field for level 0: instead write a MultiFab
        // filled with 0, with the same number of cells as the _fp field
        WriteZeroRawField( *F_fp, dm, filename, level_prefix, field_name+"_cp", lev, ng );
    } else {
        // Create an alias to the component `icomp` of F_cp
        MultiFab F_comp(*F_cp, amrex::make_alias, icomp, 1);
        // Interpolate coarse data onto fine grid
        const int r_ratio = WarpX::GetInstance().refRatio(lev-1)[0];
        const Real* dx = WarpX::GetInstance().Geom(lev-1).CellSize();
        auto F = getInterpolatedScalar( F_comp, *F_fp, dm, r_ratio, dx, ng );
        // Write interpolated raw data
        WriteRawField( *F, dm, filename, level_prefix, field_name+"_cp", lev, plot_guards );
    }
}

/** \brief Write the coarse vector multifab `F*_cp` to the file `filename`
 *  *after* sampling/interpolating its value on the fine grid corresponding
 *  to `F*_fp`. This is mainly needed because the yt reader requires the
 *  coarse and fine patch to have the same shape.
 */
void
WriteCoarseVector( const std::string field_name,
    const std::unique_ptr<MultiFab>& Fx_cp,
    const std::unique_ptr<MultiFab>& Fy_cp,
    const std::unique_ptr<MultiFab>& Fz_cp,
    const std::unique_ptr<MultiFab>& Fx_fp,
    const std::unique_ptr<MultiFab>& Fy_fp,
    const std::unique_ptr<MultiFab>& Fz_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards )
{
    int ng = 0;
    if (plot_guards) ng = Fx_fp->nGrow();

    if (lev == 0) {
        // No coarse field for level 0: instead write a MultiFab
        // filled with 0, with the same number of cells as the _fp field
        WriteZeroRawField( *Fx_fp, dm, filename, level_prefix, field_name+"x_cp", lev, ng );
        WriteZeroRawField( *Fy_fp, dm, filename, level_prefix, field_name+"y_cp", lev, ng );
        WriteZeroRawField( *Fz_fp, dm, filename, level_prefix, field_name+"z_cp", lev, ng );
    } else {
        // Interpolate coarse data onto fine grid
        const int r_ratio = WarpX::GetInstance().refRatio(lev-1)[0];
        const Real* dx = WarpX::GetInstance().Geom(lev-1).CellSize();
        auto F = getInterpolatedVector( Fx_cp, Fy_cp, Fz_cp, Fx_fp, Fy_fp, Fz_fp,
                                    dm, r_ratio, dx, ng );
        // Write interpolated raw data
        WriteRawField( *F[0], dm, filename, level_prefix, field_name+"x_cp", lev, plot_guards );
        WriteRawField( *F[1], dm, filename, level_prefix, field_name+"y_cp", lev, plot_guards );
        WriteRawField( *F[2], dm, filename, level_prefix, field_name+"z_cp", lev, plot_guards );
    }
}

/** \brief Samples/Interpolates the coarse scalar multifab `F_cp` on the
  * fine grid associated with the fine multifab `F_fp`.
  */
std::unique_ptr<MultiFab>
getInterpolatedScalar(
    const MultiFab& F_cp, const MultiFab& F_fp,
    const DistributionMapping& dm, const int r_ratio,
    const Real* /*dx*/, const int ngrow )
{
    // Prepare the structure that will contain the returned fields
    std::unique_ptr<MultiFab> interpolated_F;
    interpolated_F.reset( new MultiFab(F_fp.boxArray(), dm, 1, ngrow) );
    interpolated_F->setVal(0.);

    // Loop through the boxes and interpolate the values from the _cp data
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox ffab; // Temporary array ; contains interpolated fields
        for (MFIter mfi(*interpolated_F); mfi.isValid(); ++mfi)
        {
            Box finebx = mfi.fabbox();
            finebx.coarsen(r_ratio).refine(r_ratio); // so that finebx is coarsenable

            const FArrayBox& cfab = (F_cp)[mfi];
            ffab.resize(finebx);

            // - Fully nodal
            if ( F_fp.is_nodal() ){
                IntVect refinement_vector{AMREX_D_DECL(r_ratio, r_ratio, r_ratio)};
                node_bilinear_interp.interp(cfab, 0, ffab, 0, 1,
                        finebx, refinement_vector, {}, {}, {}, 0, 0, RunOn::Cpu);
            } else {
                amrex::Abort("Unknown field staggering.");
            }

            // Add temporary array to the returned structure
            const Box& bx = (*interpolated_F)[mfi].box();
            (*interpolated_F)[mfi].plus<RunOn::Host>(ffab, bx, bx, 0, 0, 1);
        }
    }
    return interpolated_F;
}

/** \brief Samples/Interpolates the coarse vector multifab `F*_cp` on the
  * fine grid associated with the fine multifab `F*_fp`.
  */
std::array<std::unique_ptr<MultiFab>, 3>
getInterpolatedVector(
    const std::unique_ptr<MultiFab>& Fx_cp,
    const std::unique_ptr<MultiFab>& Fy_cp,
    const std::unique_ptr<MultiFab>& Fz_cp,
    const std::unique_ptr<MultiFab>& Fx_fp,
    const std::unique_ptr<MultiFab>& Fy_fp,
    const std::unique_ptr<MultiFab>& Fz_fp,
    const DistributionMapping& dm, const int r_ratio,
    const Real* dx, const int ngrow )
{

    // Prepare the structure that will contain the returned fields
    std::array<std::unique_ptr<MultiFab>, 3> interpolated_F;
    interpolated_F[0].reset( new MultiFab(Fx_fp->boxArray(), dm, 1, ngrow) );
    interpolated_F[1].reset( new MultiFab(Fy_fp->boxArray(), dm, 1, ngrow) );
    interpolated_F[2].reset( new MultiFab(Fz_fp->boxArray(), dm, 1, ngrow) );
    for (int i=0; i<3; i++) interpolated_F[i]->setVal(0.);

    // Loop through the boxes and interpolate the values from the _cp data
    const int use_limiter = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,3> ffab; // Temporary array ; contains interpolated fields
        for (MFIter mfi(*interpolated_F[0]); mfi.isValid(); ++mfi)
        {
            Box ccbx = mfi.fabbox();
            ccbx.enclosedCells();
            ccbx.coarsen(r_ratio).refine(r_ratio); // so that ccbx is coarsenable

            const FArrayBox& cxfab = (*Fx_cp)[mfi];
            const FArrayBox& cyfab = (*Fy_cp)[mfi];
            const FArrayBox& czfab = (*Fz_cp)[mfi];
            ffab[0].resize(amrex::convert(ccbx,(*Fx_fp)[mfi].box().type()));
            ffab[1].resize(amrex::convert(ccbx,(*Fy_fp)[mfi].box().type()));
            ffab[2].resize(amrex::convert(ccbx,(*Fz_fp)[mfi].box().type()));

            // - Face centered, in the same way as B on a Yee grid
            if ( (*Fx_fp)[mfi].box().type() == IntVect{AMREX_D_DECL(1,0,0)} ){
#if (AMREX_SPACEDIM == 3)
                amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                             BL_TO_FORTRAN_ANYD(ffab[0]),
                                             BL_TO_FORTRAN_ANYD(ffab[1]),
                                             BL_TO_FORTRAN_ANYD(ffab[2]),
                                             BL_TO_FORTRAN_ANYD(cxfab),
                                             BL_TO_FORTRAN_ANYD(cyfab),
                                             BL_TO_FORTRAN_ANYD(czfab),
                                             dx, &r_ratio, &use_limiter);
#else
                amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                             BL_TO_FORTRAN_ANYD(ffab[0]),
                                             BL_TO_FORTRAN_ANYD(ffab[2]),
                                             BL_TO_FORTRAN_ANYD(cxfab),
                                             BL_TO_FORTRAN_ANYD(czfab),
                                             dx, &r_ratio, &use_limiter);
                amrex_interp_cc_bfield(ccbx.loVect(), ccbx.hiVect(),
                                       BL_TO_FORTRAN_ANYD(ffab[1]),
                                       BL_TO_FORTRAN_ANYD(cyfab),
                                       &r_ratio, &use_limiter);
#endif
            // - Edge centered, in the same way as E on a Yee grid
            } else if ( (*Fx_fp)[mfi].box().type() == IntVect{AMREX_D_DECL(0,1,1)} ){
#if (AMREX_SPACEDIM == 3)
                amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                    BL_TO_FORTRAN_ANYD(ffab[0]),
                                    BL_TO_FORTRAN_ANYD(ffab[1]),
                                    BL_TO_FORTRAN_ANYD(ffab[2]),
                                    BL_TO_FORTRAN_ANYD(cxfab),
                                    BL_TO_FORTRAN_ANYD(cyfab),
                                    BL_TO_FORTRAN_ANYD(czfab),
                                    &r_ratio, &use_limiter);
#else
                amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                    BL_TO_FORTRAN_ANYD(ffab[0]),
                                    BL_TO_FORTRAN_ANYD(ffab[2]),
                                    BL_TO_FORTRAN_ANYD(cxfab),
                                    BL_TO_FORTRAN_ANYD(czfab),
                                    &r_ratio,&use_limiter);
                amrex_interp_nd_efield(ccbx.loVect(), ccbx.hiVect(),
                                       BL_TO_FORTRAN_ANYD(ffab[1]),
                                       BL_TO_FORTRAN_ANYD(cyfab),
                                       &r_ratio);
#endif
            } else {
                amrex::Abort("Unknown field staggering.");
            }

            // Add temporary array to the returned structure
            for (int i = 0; i < 3; ++i) {
                const Box& bx = (*interpolated_F[i])[mfi].box();
                (*interpolated_F[i])[mfi].plus<RunOn::Host>(ffab[i], bx, bx, 0, 0, 1);
            }
        }
    }
    return interpolated_F;
}
