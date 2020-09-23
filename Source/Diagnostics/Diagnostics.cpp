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
#include <AMReX_Vector.H>
#include <string>
using namespace amrex::literals;

Diagnostics::Diagnostics (int i, std::string name)
    : m_diag_name(name), m_diag_index(i)
{
}

Diagnostics::~Diagnostics ()
{
    delete m_flush_format;
}

bool
Diagnostics::BaseReadParameters ()
{
    auto & warpx = WarpX::GetInstance();
    // Read list of fields requested by the user.
    amrex::ParmParse pp(m_diag_name);
    m_file_prefix = "diags/" + m_diag_name;
    pp.query("file_prefix", m_file_prefix);
    pp.query("format", m_format);
    bool varnames_specified = pp.queryarr("fields_to_plot", m_varnames);
    if (!varnames_specified){
        m_varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    }
    // If user requests rho with back-transformed diagnostics, we set plot_rho=true
    // and compute rho at each iteration
    if (WarpXUtilStr::is_in(m_varnames, "rho") && WarpX::do_back_transformed_diagnostics) {
        warpx.setplot_rho(true);
    }
    // Sanity check if user requests to plot F
    if (WarpXUtilStr::is_in(m_varnames, "F")){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            warpx.do_dive_cleaning,
            "plot F only works if warpx.do_dive_cleaning = 1");
    }
    // If user requests to plot proc_number for a serial run,
    // delete proc_number from fields_to_plot
    if (amrex::ParallelDescriptor::NProcs() == 1){
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
    // For a moving window simulation, the user-defined m_lo and m_hi must be converted.
    if (warpx.do_moving_window) {
#if (AMREX_SPACEDIM == 3)
    amrex::Vector<int> dim_map {0, 1, 2};
#else
    amrex::Vector<int> dim_map {0, 2};
#endif
       if (warpx.boost_direction[ dim_map[warpx.moving_window_dir] ] == 1) {
           // Convert user-defined lo and hi for diagnostics to account for boosted-frame
           // simulations with moving window
           amrex::Real convert_factor = 1._rt/(warpx.gamma_boost * (1._rt - warpx.beta_boost) );
           // Assuming that the window travels with speed c
           m_lo[warpx.moving_window_dir] *= convert_factor;
           m_hi[warpx.moving_window_dir] *= convert_factor;
       }
    }

    // Initialize cr_ratio with default value of 1 for each dimension.
    amrex::Vector<int> cr_ratio(AMREX_SPACEDIM, 1);
    // Read user-defined coarsening ratio for the output MultiFab.
    bool cr_specified = pp.queryarr("coarsening_ratio", cr_ratio);
    if (cr_specified) {
       for (int idim =0; idim < AMREX_SPACEDIM; ++idim) {
           m_crse_ratio[idim] = cr_ratio[idim];
       }
    }

    bool species_specified = pp.queryarr("species", m_species_names);

    // Prepare to dump rho per species:
    // - add the string "rho_<species_names>" to m_varnames if "rho" is listed in
    //   <diag_name>.<species_name>.variables in the input file
    // - while looping over all species, count the number ns_dump_rho of species
    //   that request to dump rho per species
    // - allocate and fill the array m_rho_per_species_index, which contains the
    //   indices of the species that request to dump rho per species; the indices
    //   will be then passed to the constructor of RhoFunctor and used in
    //   RhoFunctor::operator() to get the correct particle container for the given species
    amrex::Vector<std::string> species_variables;
    const MultiParticleContainer& mpc = warpx.GetPartContainer();
    // If <diag_name>.species is not provided, loop over all species to check whether
    // rho per species is requested
    if (m_species_names.size() == 0u) m_species_names = mpc.GetSpeciesNames();
    // ns_dump_rho: number of species that dump rho per species
    int ns_dump_rho = 0;
    // Loop over all species
    for (const auto& species_name : m_species_names) {
        // Parse input argument <diag_name>.<species_name>.variables
        amrex::ParmParse pp_sp(m_diag_name + "." + species_name);
        if (pp_sp.queryarr("variables", species_variables)) {
            for (const auto& var : species_variables) {
                if (var == "rho") {
                    // Add string "rho_<species_name>" to m_varnames
                    m_varnames.push_back("rho_" + species_name);
                    // Count number of species that dump rho per species
                    ns_dump_rho++;
                }
            }
        }
        // Clear content of species_variables before moving to the next species
        species_variables.clear();
    }
    // Allocate array of species indices that dump rho per species
    m_rho_per_species_index.resize(ns_dump_rho);
    // ns: total number of species
    const int ns = int(m_species_names.size());
    // is_dump_rho: species index to loop over species that dump rho per species
    int is_dump_rho = 0;
    // Loop over all species
    for (int is = 0; is < ns; is++) {
        // Parse input argument <diag_name>.<species_name>.variables
        amrex::ParmParse pp_sp(m_diag_name + "." + m_species_names[is]);
        if (pp_sp.queryarr("variables", species_variables)) {
            for (const auto& var : species_variables) {
                if (var == "rho") {
                    // Fill array of species indices that dump rho per species
                    m_rho_per_species_index.at(is_dump_rho) = is;
                    is_dump_rho++;
                }
            }
        }
        // Clear content of species_variables before moving to the next species
        species_variables.clear();
    }

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
Diagnostics::InitData ()
{
    // initialize member variables and arrays in base class::Diagnostics
    InitBaseData();
    // initialize member variables and arrays specific to each derived class
    // (FullDiagnostics, BTDiagnostics, etc.)
    DerivedInitData();
    // loop over all buffers
    for (int i_buffer = 0; i_buffer < m_num_buffers; ++i_buffer) {
        // loop over all levels
        for (int lev = 0; lev < nmax_lev; ++lev) {
            // allocate and initialize m_all_field_functors depending on diag type
            InitializeFieldFunctors(lev);
            // Initialize field buffer data, m_mf_output
            InitializeFieldBufferData(i_buffer, lev);
        }
    }
    // When particle buffers, m_particle_buffers are included, they will be initialized here
    InitializeParticleBuffer();

    amrex::ParmParse pp(m_diag_name);
    amrex::Vector <amrex::Real> dummy_val(AMREX_SPACEDIM);
    if ( pp.queryarr("diag_lo", dummy_val) || pp.queryarr("diag_hi", dummy_val) ) {
        // set geometry filter for particle-diags to true when the diagnostic domain-extent
        // is specified by the user
        for (int i = 0; i < m_all_species.size(); ++i) {
            m_all_species[i].m_do_geom_filter = true;
        }
        // Disabling particle-io for reduced domain diagnostics by reducing
        // the particle-diag vector to zero.
        // This is a temporary fix until particle_buffer is supported in diagnostics.
        m_all_species.clear();
        amrex::Print() << " WARNING: For full diagnostics on a reduced domain, particle io is not supported, yet! Therefore, particle-io is disabled for this diag " << m_diag_name << "\n";
    }
}


void
Diagnostics::InitBaseData ()
{
    auto & warpx = WarpX::GetInstance();
    // Number of levels in the simulation at the current timestep
    nlev = warpx.finestLevel() + 1;
    // default number of levels to be output = nlev
    nlev_output = nlev;
    // Maximum number of levels that will be allocated in the simulation
    nmax_lev = warpx.maxLevel() + 1;
    m_all_field_functors.resize( nmax_lev );

    // For restart, move the m_lo and m_hi of the diag consistent with the
    // current moving_window location
    if (warpx.do_moving_window) {
        int moving_dir = warpx.moving_window_dir;
        int shift_num_base = static_cast<int>((warpx.getmoving_window_x() - m_lo[moving_dir]) / warpx.Geom(0).CellSize(moving_dir) );
        m_lo[moving_dir] += shift_num_base * warpx.Geom(0).CellSize(moving_dir);
        m_hi[moving_dir] += shift_num_base * warpx.Geom(0).CellSize(moving_dir);
    }
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

    // allocate vector of geometry objects corresponding to each output multifab.
    m_geom_output.resize( m_num_buffers );
    for (int i = 0; i < m_num_buffers; ++i) {
        m_geom_output[i].resize( nmax_lev );
    }
}

void
Diagnostics::ComputeAndPack ()
{
    // prepare the field-data necessary to compute output data
    PrepareFieldDataForOutput();
    // compute the necessary fields and stiore result in m_mf_output.
    for (int i_buffer = 0; i_buffer < m_num_buffers; ++i_buffer) {
        for(int lev=0; lev<nlev_output; lev++){
            int icomp_dst = 0;
            for (int icomp=0, n=m_all_field_functors[0].size(); icomp<n; icomp++){
                // Call all functors in m_all_field_functors[lev]. Each of them computes
                // a diagnostics and writes in one or more components of the output
                // multifab m_mf_output[lev].
                m_all_field_functors[lev][icomp]->operator()(m_mf_output[i_buffer][lev], icomp_dst, i_buffer);
                // update the index of the next component to fill
                icomp_dst += m_all_field_functors[lev][icomp]->nComp();
            }
            // Check that the proper number of components of mf_avg were updated.
            AMREX_ALWAYS_ASSERT( icomp_dst == m_varnames.size() );
        }
    }
}


void
Diagnostics::FilterComputePackFlush (int step, bool force_flush)
{
    WARPX_PROFILE("Diagnostics::FilterComputePackFlush()");

    MovingWindowAndGalileanDomainShift ();

    if ( DoComputeAndPack (step, force_flush) ) {
        ComputeAndPack();

        for (int i_buffer = 0; i_buffer < m_num_buffers; ++i_buffer) {
            if ( !DoDump (step, i_buffer, force_flush) ) continue;
                Flush(i_buffer);
        }

    }

}
