#include "Diagnostics.H"
#include "WarpX.H"
#include "Average.H"
#include "Utils/WarpXUtil.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : diag_index(i), diag_name(name)
{
    ReadParameters();
}

Diagnostics::~Diagnostics ()
{
    delete m_flush_format;
}

void
Diagnostics::ReadParameters ()
{
    auto & warpx = WarpX::GetInstance();
    // Read list of fields requested by the user.
    ParmParse pp(diag_name);
    pp.query("file_prefix", file_prefix);
    pp.query("period", m_period);
    pp.query("plot_raw_fields", m_plot_raw_fields);
    pp.query("plot_raw_fields_guards", m_plot_raw_fields_guards);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_plot_F==false, "cannot plot_F yet");
    if (!pp.queryarr("fields_to_plot", varnames)){
        varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    }
    ncomp = varnames.size();
    // set plot_rho to true of the users requests it, so that
    // rho is computed at each iteration.
    if (WarpXUtilStr::is_in(varnames, "rho")) warpx.setplot_rho(true);
    // Sanity check if user requests to plot F
    if (WarpXUtilStr::is_in(varnames, "F")){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            warpx.do_dive_cleaning,
            "plot F only works if warpx.do_dive_cleaning = 1");
    }
    // If user requests to plot proc_number for a serial run,
    // delete proc_number from fields_to_plot
    if (ParallelDescriptor::NProcs() == 1){
        varnames.erase(
            std::remove(varnames.begin(), varnames.end(), "proc_number"),
            varnames.end());
    }
}

void
Diagnostics::InitData ()
{
    Print()<<"Diagnostics::InitData\n";
    auto & warpx = WarpX::GetInstance();
    nlev = warpx.finestLevel() + 1;
    // Initialize vector of pointers to the fields requested by the user.
    all_field_functors.resize( nlev );
    mf_avg.resize( nlev );
    for ( int lev=0; lev<nlev; lev++ ){
        all_field_functors[lev].resize( ncomp );
        for (int comp=0; comp<ncomp; comp++){
            if        ( varnames[comp] == "Ex" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 0));
            } else if ( varnames[comp] == "Ey" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 1));
            } else if ( varnames[comp] == "Ez" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 2));
            } else if ( varnames[comp] == "Bx" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 0));
            } else if ( varnames[comp] == "By" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 1));
            } else if ( varnames[comp] == "Bz" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 2));
            } else if ( varnames[comp] == "jx" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 0));
            } else if ( varnames[comp] == "jy" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 1));
            } else if ( varnames[comp] == "jz" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 2));
            } else if ( varnames[comp] == "rho" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_rho_fp(lev));
            } else if ( varnames[comp] == "part_per_cell" ){
                amrex::Abort("plotting part_per_cell is not supported");
            } else {
                amrex::Abort("Unknown field in fields_to_plot");
            }
        }
        // Allocate output multifab
        // Note: default MultiFab constructor is cell-centered
        mf_avg[lev] = MultiFab(warpx.boxArray(lev),
                               warpx.DistributionMap(lev),
                               ncomp, 0);
    }
    // Construct Flush class. So far, only Plotfile is implemented.
    m_flush_format = new FlushFormatPlotfile;
}

void
Diagnostics::ComputeAndPack ()
{
    // First, make sure all guard cells are properly filled
    // Probably overkill/unnecessary, but safe and shouldn't happen often !!
    auto & warpx = WarpX::GetInstance();
    warpx.FillBoundaryE(warpx.getngE(), warpx.getngExtra());
    warpx.FillBoundaryB(warpx.getngE(), warpx.getngExtra());
#ifndef WARPX_USE_PSATD
    warpx.FillBoundaryAux(warpx.getngUpdateAux());
#endif
    warpx.UpdateAuxilaryData();

    warpx.FieldGather();

    // cell-center fields and store result in mf_avg.
    for(int lev=0; lev<nlev; lev++){
        for (int icomp=0; icomp<ncomp; icomp++){
            all_field_functors[lev][icomp]->operator()(mf_avg[lev], icomp, 1);
        }
    }
}

void
Diagnostics::Flush ()
{
    auto & warpx = WarpX::GetInstance();
    m_flush_format->WriteToFile(
        varnames, GetVecOfConstPtrs(mf_avg), warpx.Geom(), warpx.getistep(),
        warpx.gett_new(0), warpx.GetPartContainer(), nlev, file_prefix,
        m_plot_raw_fields, m_plot_raw_fields_guards, m_plot_rho, m_plot_F);
}

void
Diagnostics::FlushRaw () {}

bool
Diagnostics::DoDump (int step, bool force_flush)
{
    if (force_flush) return true;
    if ( m_period>0 && (step+1)%m_period==0 ) return true;
    return false;
}

void
Diagnostics::AddToVarNames (std::string name, std::string suffix) {
    auto coords = {"x", "y", "z"};
    for(auto coord:coords) varnames.push_back(name+coord+suffix);
    ncomp += 1;
}
