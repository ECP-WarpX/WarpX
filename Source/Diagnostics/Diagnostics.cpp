
#include "Diagnostics.H"
#include "WarpX.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : diag_index(i), diag_name(name)
{
}

void
Diagnostics::InitData ()
{
    Print()<<"Diagnostics::InitData\n";

    auto & warpx = WarpX::GetInstance();
    finest_level = warpx.finestLevel();
    
    allfields.resize( finest_level );
    for ( int lev=0; lev<finest_level; lev++ ){
        allfields[lev].resize( ncomp );
        for ( int dim=0; dim<3; dim++ ){
            allfields[lev][dim  ] = warpx.get_pointer_Efield_aux(lev, dim);
            allfields[lev][dim+3] = warpx.get_pointer_Bfield_aux(lev, dim);
            allfields[lev][dim+6] = warpx.get_pointer_current_fp(lev, dim);
        }
    }
}

void
Diagnostics::Filter () {}

void
Diagnostics::PackFields ()
{
    // Average the fields from the simulation grid to the cell centers
    const int ngrow = 0;
    
    AverageAndPackFields( varnames, mf_avg, ngrow );
    output_mf = amrex::GetVecOfConstPtrs(mf_avg);
    output_geom = Geom();
}

void
Diagnostics::Flush () {}

void
Diagnostics::FlushRaw () {}

void
Diagnostics::AverageAndPackFields (const int ngrow) const
{
    const int nvecs = 3;
    auto & warpx = WarpX::GetInstance();
    // Loop over levels of refinement
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Allocate pointers to the `ncomp` fields that will be added
        mf_avg.push_back( MultiFab(warpx.grids[lev], warpx.dmap[lev], ncomp, ngrow));

        // For E, B and J, if at least one component is requested,
        // build cell-centered temporary MultiFab with 3 comps
        MultiFab mf_tmp_E, mf_tmp_B, mf_tmp_J;
        // Build mf_tmp_E is at least one component of E is requested
        // Allocate temp MultiFab with 3 components
        mf_tmp_E = MultiFab(warpx.grids[lev], warpx.dmap[lev], nvecs, ngrow);
        // Fill MultiFab mf_tmp_E with averaged E
        AverageAndPackVectorField(mf_tmp_E, Efield_aux[lev], dmap[lev], 0, ngrow);
        int dcomp = 3;
        AverageAndPackVectorFieldComponents(mf_tmp_E, Efield_aux[lev], dmap[lev], dcomp, ngrow);
        // Same for B
        mf_tmp_B = MultiFab(grids[lev], dmap[lev], nvecs, ngrow);
        AverageAndPackVectorField(mf_tmp_B, Bfield_aux[lev], dmap[lev], 0, ngrow);
        dcomp = 3;
        AverageAndPackVectorFieldComponents(mf_tmp_B, Bfield_aux[lev], dmap[lev], dcomp, ngrow);
        // Same for J
        mf_tmp_J = MultiFab(grids[lev], dmap[lev], nvecs, ngrow);
        AverageAndPackVectorField(mf_tmp_J, current_fp[lev], dmap[lev], 0, ngrow);
        dcomp = 3;
        AverageAndPackVectorFieldComponents(mf_tmp_J, current_fp[lev], dmap[lev], dcomp, ngrow);
        
        int dcomp;
        // Go through the different fields in fields_to_plot, pack them into
        // mf_avg[lev] add the corresponding names to `varnames`.
        // plot_fine_patch and plot_coarse_patch are treated separately
        // (after this for loop).
        dcomp = 0;
        for (int ifield=0; ifield<ncomp; ifield++){
            std::string fieldname = varnames[ifield];
            if        (fieldname == "Ex"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_E, 0, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_E, 0, dcomp, ngrow, "Er", varnames);
            } else if (fieldname == "Ey"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_E, 1, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_E, 1, dcomp, ngrow, "Etheta", varnames);
            } else if (fieldname == "Ez"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_E, 2, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_E, 2, dcomp, ngrow, "Ez", varnames);
            } else if (fieldname == "Bx"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_B, 0, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_B, 0, dcomp, ngrow, "Br", varnames);
            } else if (fieldname == "By"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_B, 1, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_B, 1, dcomp, ngrow, "Btheta", varnames);
            } else if (fieldname == "Bz"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_B, 2, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_B, 2, dcomp, ngrow, "Bz", varnames);
            } else if (fieldname == "jx"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_J, 0, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_J, 0, dcomp, ngrow, "jr", varnames);
            } else if (fieldname == "jy"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_J, 1, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_J, 1, dcomp, ngrow, "jtheta", varnames);
            } else if (fieldname == "jz"){
                MultiFab::Copy( mf_avg[lev], mf_tmp_J, 2, dcomp++, 1, ngrow);
                CopyVectorFieldComponentsToMultiFab(lev, mf_avg, mf_tmp_J, 2, dcomp, ngrow, "jz", varnames);
            } else {
                amrex::Abort("unknown field in fields_to_plot: " + fieldname);
            }
        }
        BL_ASSERT(dcomp == ncomp);
    } // end loop over levels of refinement
};

//void
//AverageAndPackVectorField( MultiFab& mf_avg,
//                           const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field,
//                           const DistributionMapping& dm,
//                           const int dcomp, const int ngrow )
void
AverageAndPackVectorField( const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field,
                           const DistributionMapping& dm,
                           const int dcomp, const int ngrow )
{
#ifndef WARPX_DIM_RZ
    (void)dm;
#endif
    // The object below is temporary, and is needed because
    // `average_edge_to_cellcenter` requires fields to be passed as Vector
    Vector<const MultiFab*> srcmf(AMREX_SPACEDIM);

    // Check the type of staggering of the 3-component `vector_field`
    // and average accordingly:
    // - Fully cell-centered field (no average needed; simply copy)
    if ( vector_field[0]->is_cell_centered() ){

        MultiFab::Copy( mf_avg, *vector_field[0], 0, dcomp  , 1, ngrow);
        MultiFab::Copy( mf_avg, *vector_field[1], 0, dcomp+1, 1, ngrow);
        MultiFab::Copy( mf_avg, *vector_field[2], 0, dcomp+2, 1, ngrow);

        // - Fully nodal
    } else if ( vector_field[0]->is_nodal() ){

        amrex::average_node_to_cellcenter( mf_avg, dcomp  ,
                                          *vector_field[0], 0, 1, ngrow);
        amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                          *vector_field[1], 0, 1, ngrow);
        amrex::average_node_to_cellcenter( mf_avg, dcomp+2,
                                          *vector_field[2], 0, 1, ngrow);

        // - Face centered, in the same way as B on a Yee grid
    } else if ( vector_field[0]->is_nodal(0) ){

        // Note that average_face_to_cellcenter operates only on the number of
        // arrays equal to the number of dimensions. So, for 2D, PackPlotDataPtrs
        // packs in the x and z (or r and z) arrays, which are then cell averaged.
        // The Copy code then copies the z from the 2nd to the 3rd field,
        // and copies over directly the y (or theta) component (which is
        // already cell centered).
        if (vector_field[0]->nComp() > 1) {
#ifdef WARPX_DIM_RZ
            // When there are more than one components, the total
            // fields needs to be constructed in temporary MultiFabs.
            // Note that mf_total is declared in the same way as
            // vector_field so that it can be passed into PackPlotDataPtrs.
            std::array<std::unique_ptr<MultiFab>,3> mf_total;
            mf_total[0].reset(new MultiFab(vector_field[0]->boxArray(), dm, 1, vector_field[0]->nGrowVect()));
            mf_total[1].reset(new MultiFab(vector_field[1]->boxArray(), dm, 1, vector_field[1]->nGrowVect()));
            mf_total[2].reset(new MultiFab(vector_field[2]->boxArray(), dm, 1, vector_field[2]->nGrowVect()));
            ConstructTotalRZField(mf_total, vector_field);
            PackPlotDataPtrs(srcmf, mf_total);
            amrex::average_face_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            MultiFab::Copy( mf_avg, *mf_total[1], 0, dcomp+1, 1, ngrow);
#else
           amrex::Abort("AverageAndPackVectorField not implemented for ncomp > 1");
#endif
        } else {
            PackPlotDataPtrs(srcmf, vector_field);
            amrex::average_face_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
#if (AMREX_SPACEDIM == 2)
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            MultiFab::Copy( mf_avg, *vector_field[1], 0, dcomp+1, 1, ngrow);
#endif
        }

        // - Edge centered, in the same way as E on a Yee grid
    } else if ( !vector_field[0]->is_nodal(0) ){

        // See comment above, though here, the y (or theta) component
        // has node centering.
        if (vector_field[0]->nComp() > 1) {
#ifdef WARPX_DIM_RZ
            // When there are more than one components, the total
            // fields needs to be constructed in temporary MultiFabs
            // Note that mf_total is declared in the same way as
            // vector_field so that it can be passed into PackPlotDataPtrs.
            std::array<std::unique_ptr<MultiFab>,3> mf_total;
            mf_total[0].reset(new MultiFab(vector_field[0]->boxArray(), dm, 1, vector_field[0]->nGrowVect()));
            mf_total[1].reset(new MultiFab(vector_field[1]->boxArray(), dm, 1, vector_field[1]->nGrowVect()));
            mf_total[2].reset(new MultiFab(vector_field[2]->boxArray(), dm, 1, vector_field[2]->nGrowVect()));
            ConstructTotalRZField(mf_total, vector_field);
            PackPlotDataPtrs(srcmf, mf_total);
            amrex::average_edge_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                              *mf_total[1], 0, 1, ngrow);
#else
           amrex::Abort("AverageAndPackVectorField not implemented for ncomp > 1");
#endif
        } else {
            PackPlotDataPtrs(srcmf, vector_field);
            amrex::average_edge_to_cellcenter( mf_avg, dcomp, srcmf, ngrow);
#if (AMREX_SPACEDIM == 2)
            MultiFab::Copy( mf_avg, mf_avg, dcomp+1, dcomp+2, 1, ngrow);
            amrex::average_node_to_cellcenter( mf_avg, dcomp+1,
                                              *vector_field[1], 0, 1, ngrow);
#endif
        }

    } else {
        amrex::Abort("Unknown staggering.");
    }
}
