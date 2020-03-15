
#include "Diagnostics.H"
#include "WarpX.H"
#include "Average.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : diag_index(i), diag_name(name)
{}

void
Diagnostics::InitData ()
{
    Print()<<"Diagnostics::InitData\n";
    
    auto & warpx = WarpX::GetInstance();
    nlev = warpx.finestLevel() + 1;
    allfields.resize( nlev );
    mf_avg.resize( nlev );
    for ( int lev=0; lev<nlev; lev++ ){
        allfields[lev].resize( ncomp );
        for ( int dim=0; dim<3; dim++ ){
            allfields[lev][dim  ] = warpx.get_pointer_Efield_aux(lev, dim);
            allfields[lev][dim+3] = warpx.get_pointer_Bfield_aux(lev, dim);
            allfields[lev][dim+6] = warpx.get_pointer_current_fp(lev, dim);
        }
        // Default MultiFab constructor -> cell-centered
        mf_avg[lev] = MultiFab(warpx.boxArray(lev),
                               warpx.DistributionMap(lev),
                               ncomp, 0);
    }
}

void
Diagnostics::FilterAndPack ()
{
    for(int lev=0; lev<nlev; lev++){
        for (int icomp=0; icomp<ncomp; icomp++){
            Average::ToCellCenter ( mf_avg[lev],
                                    *allfields[lev][icomp],
                                    icomp, 0 );
        }
    }
}

void
Diagnostics::Flush () {
    auto & warpx = WarpX::GetInstance();
    //const auto step = istep[0];
    //const std::string& plotfilename = amrex::Concatenate(plot_file,step);
    const std::string& plotfilename = "toto";
    amrex::Print() << "  Writing plotfile " << plotfilename << "\n";

    //Vector<std::string> varnames; // Name of the written fields
    //Vector<MultiFab> mf_avg; // contains the averaged, cell-centered fields
    //Vector<const MultiFab*> output_mf; // will point to the data to be written
    //Vector<Geometry> output_geom;

    //prepareFields(step, varnames, mf_avg, output_mf, output_geom);

    // Write the fields contained in `mf_avg`, and corresponding to the
    // names `varnames`, into a plotfile.
    // Prepare extra directory (filled later), for the raw fields
    Vector<std::string> rfs;
    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::Version_v1);
    amrex::WriteMultiLevelPlotfile(plotfilename, nlev,
                                   GetVecOfConstPtrs(mf_avg),
                                   varnames, warpx.Geom(),
                                   0., warpx.getistep(), warpx.refRatio(),
                                   "HyperCLaw-V1.1",
                                   "Level_",
                                   "Cell",
                                   rfs
                                   );

    // mypc->WritePlotFile(plotfilename);

    // WriteJobInfo(plotfilename);

    // WriteWarpXHeader(plotfilename);

    // VisMF::SetHeaderVersion(current_version);
}

void
Diagnostics::FlushRaw () {}
