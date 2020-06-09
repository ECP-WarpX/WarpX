#include "Diagnostics.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/PartPerCellFunctor.H"
#include "ComputeDiagFunctors/PartPerGridFunctor.H"
#include "ComputeDiagFunctors/DivBFunctor.H"
#include "ComputeDiagFunctors/DivEFunctor.H"
#include "FlushFormats/FlushFormatPlotfile.H"
#include "FlushFormats/FlushFormatCheckpoint.H"
#include "FlushFormats/FlushFormatAscent.H"
#include "FlushFormats/FlushFormatSensei.H"
#ifdef WARPX_USE_OPENPMD
#   include "FlushFormats/FlushFormatOpenPMD.H"
#endif
#include "WarpX.H"
#include "Utils/WarpXUtil.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : m_diag_name(name), m_diag_index(i)
{
}

Diagnostics::~Diagnostics ()
{
    delete m_flush_format;
}

bool
Diagnostics::ReadBaseParameters ()
{
    auto & warpx = WarpX::GetInstance();
    // Read list of fields requested by the user.
    ParmParse pp(m_diag_name);
    m_file_prefix = "diags/" + m_diag_name;
    pp.query("file_prefix", m_file_prefix);
    pp.query("format", m_format);
    bool varnames_specified = pp.queryarr("fields_to_plot", m_varnames);
    if (!varnames_specified){
        m_varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    }
    // set plot_rho to true of the users requests it, so that
    // rho is computed at each iteration.
    if (WarpXUtilStr::is_in(m_varnames, "rho")) warpx.setplot_rho(true);
    // Sanity check if user requests to plot F
    if (WarpXUtilStr::is_in(m_varnames, "F")){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            warpx.do_dive_cleaning,
            "plot F only works if warpx.do_dive_cleaning = 1");
    }
    // If user requests to plot proc_number for a serial run,
    // delete proc_number from fields_to_plot
    if (ParallelDescriptor::NProcs() == 1){
        m_varnames.erase(
            std::remove(m_varnames.begin(), m_varnames.end(), "proc_number"),
            m_varnames.end());
    }

    // Read user-defined physical extents for the output and store in m_lo and m_hi.
    m_lo.resize(AMREX_SPACEDIM);
    m_hi.resize(AMREX_SPACEDIM);

    bool lo_specified = pp.queryarr("diag_lo", m_lo);

    if (!lo_specified) {
       for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
            m_lo[idim] = warpx.Geom(0).ProbLo(idim);
       }
    }
    bool hi_specified = pp.queryarr("diag_hi", m_hi);
    if (!hi_specified) {
       for (int idim =0; idim < AMREX_SPACEDIM; ++idim) {
            m_hi[idim] = warpx.Geom(0).ProbHi(idim);
       }
    }

    // Initialize cr_ratio with default value of 1 for each dimension.
    Vector<int> cr_ratio(AMREX_SPACEDIM, 1);
    // Read user-defined coarsening ratio for the output MultiFab.
    bool cr_specified = pp.queryarr("coarsening_ratio", cr_ratio);
    if (cr_specified) {
       for (int idim =0; idim < AMREX_SPACEDIM; ++idim) {
           m_crse_ratio[idim] = cr_ratio[idim];
       }
    }

    bool species_specified = pp.queryarr("species", m_species_names);

    bool checkpoint_compatibility = false;
    if (m_format == "checkpoint"){
       if ( varnames_specified == false &&
            lo_specified == false &&
            hi_specified == false &&
            cr_specified == false &&
            species_specified == false ) checkpoint_compatibility = true;
    }
    return checkpoint_compatibility;

}

void
Diagnostics::InitBaseData ()
{
    auto & warpx = WarpX::GetInstance();
    // Number of levels
    nlev = warpx.finestLevel() + 1;
    // Maximum number of levels that will be allocated in the simulation
    nmax_lev = warpx.maxLevel() + 1;
    m_all_field_functors.resize( nmax_lev );

    // Construct Flush class.
    if        (m_format == "plotfile"){
        m_flush_format = new FlushFormatPlotfile;
    } else if (m_format == "checkpoint"){
        // creating checkpoint format 
        m_flush_format = new FlushFormatCheckpoint;
    } else if (m_format == "ascent"){
        m_flush_format = new FlushFormatAscent;
    } else if (m_format == "sensei"){
#ifdef BL_USE_SENSEI_INSITU
        m_flush_format = new FlushFormatSensei(
            dynamic_cast<amrex::AmrMesh*>(const_cast<WarpX*>(&warpx)),
            m_diag_name);
#else
        amrex::Abort("To use SENSEI in situ, compile with USE_SENSEI=TRUE");
#endif
    } else if (m_format == "openpmd"){
#ifdef WARPX_USE_OPENPMD
        m_flush_format = new FlushFormatOpenPMD(m_diag_name);
#else
        amrex::Abort("To use openpmd output format, need to compile with USE_OPENPMD=TRUE");
#endif
    } else {
        amrex::Abort("unknown output format");
    }

    // allocate vector of buffers then allocate vector of levels for each buffer
    m_mf_output.resize( m_num_buffers );
    for (int i = 0; i < m_num_buffers; ++i) {
        m_mf_output[i].resize( nmax_lev );
    }
}

