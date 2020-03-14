#include "MultiDiagnostics.H"
#include <AMReX_ParmParse.H>

using namespace amrex;

MultiDiagnostics::MultiDiagnostics ()
{
    ReadParameters();
    Print()<<"ndiags"<<ndiags<<'\n';
    alldiags.resize( ndiags );
    for (int i=0; i<ndiags; i++){
        if ( diags_types[i] == DiagTypes::Full ){
            alldiags[i].reset( new FullDiagnostics(i, diags_names[i]) );
        } else {
            amrex::Abort("Unknown diagnostic type");
        }
    }
}

void
MultiDiagnostics::InitData ()
{
    for( auto& diag : alldiags ){
        diag->InitData();
    }
}

void
MultiDiagnostics::ReadParameters ()
{
    ParmParse pp("diagnostics");
    pp.query("ndiags", ndiags);
    diags_types.resize( ndiags );
    if (ndiags > 0) pp.getarr("diags_names", diags_names);
    for (int i=0; i<ndiags; i++){
        ParmParse ppd(diags_names[i]);
        std::string diag_type_str;
        ppd.get("diag_type", diag_type_str);
        Print()<<diag_type_str<<'\n';
        if (diag_type_str == "Full") diags_types[i] = DiagTypes::Full;
    }
}

void
MultiDiagnostics::FilterComputePackFlush ()
{
    for (auto& diag : alldiags){
        diag->FilterAndPack();
        diag->Flush();
        diag->FlushRaw();
    }
}
