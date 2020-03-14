
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

    }
}

void
Diagnostics::Flush () {}

void
Diagnostics::FlushRaw () {}
