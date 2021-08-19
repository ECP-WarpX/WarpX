#include "FullDiagnostics.H"

#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/DivBFunctor.H"
#include "ComputeDiagFunctors/DivEFunctor.H"
#include "ComputeDiagFunctors/PartPerCellFunctor.H"
#include "ComputeDiagFunctors/PartPerGridFunctor.H"
#include "ComputeDiagFunctors/RhoFunctor.H"
#include "Diagnostics/Diagnostics.H"
#include "Diagnostics/ParticleDiag/ParticleDiag.H"
#include "FlushFormats/FlushFormat.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_CoordSys.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MakeType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace amrex::literals;

FullDiagnostics::FullDiagnostics (int i, std::string name)
    : Diagnostics(i, name)
{
    ReadParameters();
    BackwardCompatibility();
}

void
FullDiagnostics::InitializeParticleBuffer ()
{
    // When particle buffers are included, the vector of particle containers
    // must be allocated in this function.
    // Initialize data in the base class Diagnostics
    auto & warpx = WarpX::GetInstance();

    const MultiParticleContainer& mpc = warpx.GetPartContainer();
    // If not specified, dump all species
    if (m_output_species_names.empty()) m_output_species_names = mpc.GetSpeciesNames();
    // Initialize one ParticleDiag per species requested
    for (auto const& species : m_output_species_names){
        const int idx = mpc.getSpeciesID(species);
        m_output_species.push_back(ParticleDiag(m_diag_name, species, mpc.GetParticleContainerPtr(idx)));
    }
}

void
FullDiagnostics::ReadParameters ()
{
    // Read list of full diagnostics fields requested by the user.
    bool checkpoint_compatibility = BaseReadParameters();
    amrex::ParmParse pp_diag_name(m_diag_name);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "plotfile" || m_format == "openpmd" ||
        m_format == "checkpoint" || m_format == "ascent" ||
        m_format == "sensei",
        "<diag>.format must be plotfile or openpmd or checkpoint or ascent or sensei");
    std::vector<std::string> intervals_string_vec = {"0"};
    pp_diag_name.queryarr("intervals", intervals_string_vec);
    m_intervals = IntervalsParser(intervals_string_vec);
    bool raw_specified = pp_diag_name.query("plot_raw_fields", m_plot_raw_fields);
    raw_specified += pp_diag_name.query("plot_raw_fields_guards", m_plot_raw_fields_guards);
    raw_specified += pp_diag_name.query("plot_raw_rho", m_plot_raw_rho);

#ifdef WARPX_DIM_RZ
    pp_diag_name.query("dump_rz_modes", m_dump_rz_modes);
#else
    amrex::ignore_unused(m_dump_rz_modes);
#endif

    if (m_format == "checkpoint"){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            raw_specified == false &&
            checkpoint_compatibility == true,
            "For a checkpoint output, cannot specify these parameters as all data must be dumped "
            "to file for a restart");
    }
    // Number of buffers = 1 for FullDiagnostics.
    // It is used to allocate the number of output multi-level MultiFab, m_mf_output
    m_num_buffers = 1;
}

void
FullDiagnostics::BackwardCompatibility ()
{
    amrex::ParmParse pp_diag_name(m_diag_name);
    std::vector<std::string> backward_strings;
    if (pp_diag_name.queryarr("period", backward_strings)){
        amrex::Abort("<diag_name>.period is no longer a valid option. "
                     "Please use the renamed option <diag_name>.intervals instead.");
    }
}

void
FullDiagnostics::Flush ( int i_buffer )
{
    // This function should be moved to Diagnostics when plotfiles/openpmd format
    // is supported for BackTransformed Diagnostics, in BTDiagnostics class.
    auto & warpx = WarpX::GetInstance();

    m_flush_format->WriteToFile(
        m_varnames, m_mf_output[i_buffer], m_geom_output[i_buffer], warpx.getistep(),
        warpx.gett_new(0), m_output_species, nlev_output, m_file_prefix, m_file_min_digits,
        m_plot_raw_fields, m_plot_raw_fields_guards, m_plot_raw_rho, m_plot_raw_F);

    FlushRaw();
}

void
FullDiagnostics::FlushRaw () {}


bool
FullDiagnostics::DoDump (int step, int /*i_buffer*/, bool force_flush)
{
    if (m_already_done) return false;
    if ( force_flush || (m_intervals.contains(step+1)) ){
        m_already_done = true;
        return true;
    }
    return false;
}

bool
FullDiagnostics::DoComputeAndPack (int step, bool force_flush)
{
    // Data must be computed and packed for full diagnostics
    // whenever the data needs to be flushed.
    if (force_flush || m_intervals.contains(step+1) ){
        return true;
    }
    return false;
}


void
FullDiagnostics::AddRZModesToDiags (int lev)
{
#ifdef WARPX_DIM_RZ

    if (!m_dump_rz_modes) return;

    auto & warpx = WarpX::GetInstance();
    int ncomp_multimodefab = warpx.get_pointer_Efield_aux(0, 0)->nComp();
    // Make sure all multifabs have the same number of components
    for (int dim=0; dim<3; dim++){
        AMREX_ALWAYS_ASSERT(
            warpx.get_pointer_Efield_aux(lev, dim)->nComp() == ncomp_multimodefab );
        AMREX_ALWAYS_ASSERT(
            warpx.get_pointer_Bfield_aux(lev, dim)->nComp() == ncomp_multimodefab );
        AMREX_ALWAYS_ASSERT(
            warpx.get_pointer_current_fp(lev, dim)->nComp() == ncomp_multimodefab );
    }

    // Check if divE is requested
    // If so, all components will be written out
    bool divE_requested = false;
    for (int comp = 0; comp < m_varnames.size(); comp++) {
        if ( m_varnames[comp] == "divE" ) {
            divE_requested = true;
        }
    }

    // If rho is requested, all components will be written out
    bool rho_requested = WarpXUtilStr::is_in( m_varnames, "rho" );

    // First index of m_all_field_functors[lev] where RZ modes are stored
    int icomp = m_all_field_functors[0].size();
    const std::array<std::string, 3> coord {"r", "theta", "z"};

    // Er, Etheta, Ez, Br, Btheta, Bz, jr, jtheta, jz
    // Each of them being a multi-component multifab
    int n_new_fields = 9;
    if (divE_requested) {
        n_new_fields += 1;
    }
    if (rho_requested) {
        n_new_fields += 1;
    }
    m_all_field_functors[lev].resize( m_all_field_functors[0].size() + n_new_fields );
    // E
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        m_all_field_functors[lev][icomp] =
            std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        AddRZModesToOutputNames(std::string("E") + coord[dim],
                                warpx.get_pointer_Efield_aux(0, 0)->nComp());
        icomp += 1;
    }
    // B
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        m_all_field_functors[lev][icomp] =
            std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        AddRZModesToOutputNames(std::string("B") + coord[dim],
                                warpx.get_pointer_Bfield_aux(0, 0)->nComp());
        icomp += 1;
    }
    // j
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        m_all_field_functors[lev][icomp] =
            std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        icomp += 1;
        AddRZModesToOutputNames(std::string("J") + coord[dim],
                                warpx.get_pointer_current_fp(0, 0)->nComp());
    }
    // divE
    if (divE_requested) {
        m_all_field_functors[lev][icomp] = std::make_unique<DivEFunctor>(warpx.get_array_Efield_aux(lev), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        icomp += 1;
        AddRZModesToOutputNames(std::string("divE"), ncomp_multimodefab);
    }
    // rho
    if (rho_requested) {
        m_all_field_functors[lev][icomp] = std::make_unique<RhoFunctor>(lev, m_crse_ratio, false, ncomp_multimodefab);
        icomp += 1;
        AddRZModesToOutputNames(std::string("rho"), ncomp_multimodefab);
    }
    // Sum the number of components in input vector m_all_field_functors
    // and check that it corresponds to the number of components in m_varnames
    // and m_mf_output
    int ncomp_from_src = 0;
    for (int i=0; i<m_all_field_functors[0].size(); i++){
        ncomp_from_src += m_all_field_functors[lev][i]->nComp();
    }
    AMREX_ALWAYS_ASSERT( ncomp_from_src == m_varnames.size() );
#else
    amrex::ignore_unused(lev);
#endif
}

void
FullDiagnostics::AddRZModesToOutputNames (const std::string& field, int ncomp){
#ifdef WARPX_DIM_RZ
    // In cylindrical geometry, real and imag part of each mode are also
    // dumped to file separately, so they need to be added to m_varnames
    m_varnames.push_back( field + "_0_real" );
    for (int ic=1 ; ic < (ncomp+1)/2 ; ic += 1) {
        m_varnames.push_back( field + "_" + std::to_string(ic) + "_real" );
        m_varnames.push_back( field + "_" + std::to_string(ic) + "_imag" );
    }
#else
    amrex::ignore_unused(field, ncomp);
#endif
}


void
FullDiagnostics::InitializeFieldBufferData (int i_buffer, int lev ) {
    auto & warpx = WarpX::GetInstance();
    amrex::RealBox diag_dom;
    bool use_warpxba = true;
    const amrex::IntVect blockingFactor = warpx.blockingFactor( lev );

    // Default BoxArray and DistributionMap for initializing the output MultiFab, m_mf_output.
    amrex::BoxArray ba = warpx.boxArray(lev);
    amrex::DistributionMapping dmap = warpx.DistributionMap(lev);
    // Check if warpx BoxArray is coarsenable.
    if (warpx.get_numprocs() == 0)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE (
            ba.coarsenable(m_crse_ratio), "Invalid coarsening ratio for field diagnostics."
            "Must be an integer divisor of the blocking factor.");
    }
    else
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE (
            ba.coarsenable(m_crse_ratio), "Invalid coarsening ratio for field diagnostics."
            "The total number of cells must be a multiple of the coarsening ratio multiplied by numprocs.");
    }

    // Find if user-defined physical dimensions are different from the simulation domain.
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
         // To ensure that the diagnostic lo and hi are within the domain defined at level, lev.
        diag_dom.setLo(idim, std::max(m_lo[idim],warpx.Geom(lev).ProbLo(idim)) );
        diag_dom.setHi(idim, std::min(m_hi[idim],warpx.Geom(lev).ProbHi(idim)) );
        if ( fabs(warpx.Geom(lev).ProbLo(idim) - diag_dom.lo(idim))
                               >  warpx.Geom(lev).CellSize(idim) )
             use_warpxba = false;
        if ( fabs(warpx.Geom(lev).ProbHi(idim) - diag_dom.hi(idim))
                               > warpx.Geom(lev).CellSize(idim) )
             use_warpxba = false;

        // User-defined value for coarsening should be an integer divisor of
        // blocking factor at level, lev. This assert is not relevant and thus
        // removed if warpx.numprocs is used for the domain decomposition.
        if (warpx.get_numprocs() == 0)
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( blockingFactor[idim] % m_crse_ratio[idim]==0,
                           " coarsening ratio must be integer divisor of blocking factor");
        }
    }

    if (use_warpxba == false) {
        // Following are the steps to compute the lo and hi index corresponding to user-defined
        // m_lo and m_hi using the same resolution as the simulation at level, lev.
        amrex::IntVect lo(0);
        amrex::IntVect hi(1);
        for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
            // lo index with same cell-size as simulation at level, lev.
            lo[idim] = std::max( static_cast<int>( floor (
                          ( diag_dom.lo(idim) - warpx.Geom(lev).ProbLo(idim)) /
                            warpx.Geom(lev).CellSize(idim)) ), 0 );
            // hi index with same cell-size as simulation at level, lev.
            hi[idim] = std::max( static_cast<int> ( ceil (
                          ( diag_dom.hi(idim) - warpx.Geom(lev).ProbLo(idim)) /
                            warpx.Geom(lev).CellSize(idim) ) ), 0) - 1 ;
            // if hi<=lo, then hi = lo + 1, to ensure one cell in that dimension
            if ( hi[idim] <= lo[idim] ) {
                 hi[idim]  = lo[idim] + 1;
                 AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    m_crse_ratio[idim]==1, "coarsening ratio in reduced dimension must be 1."
                 );
            }
        }

        // Box for the output MultiFab corresponding to the user-defined physical co-ordinates at lev.
        amrex::Box diag_box( lo, hi );
        // Define box array
        amrex::BoxArray diag_ba;
        diag_ba.define(diag_box);
        ba = diag_ba.maxSize( warpx.maxGridSize( lev ) );
        // At this point in the code, the BoxArray, ba, is defined with the same index space and
        // resolution as the simulation, at level, lev.
        // Coarsen and refine so that the new BoxArray is coarsenable.
        ba.coarsen(m_crse_ratio).refine(m_crse_ratio);

        // Update the physical co-ordinates m_lo and m_hi using the final index values
        // from the coarsenable, cell-centered BoxArray, ba.
        for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            diag_dom.setLo( idim, warpx.Geom(lev).ProbLo(idim) +
                ba.getCellCenteredBox(0).smallEnd(idim) * warpx.Geom(lev).CellSize(idim));
            diag_dom.setHi( idim, warpx.Geom(lev).ProbLo(idim) +
                (ba.getCellCenteredBox( ba.size()-1 ).bigEnd(idim) + 1) * warpx.Geom(lev).CellSize(idim));
        }
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_crse_ratio.min() > 0, "Coarsening ratio must be non-zero.");
    // The BoxArray is coarsened based on the user-defined coarsening ratio.
    ba.coarsen(m_crse_ratio);
    // Generate a new distribution map if the physical m_lo and m_hi for the output
    // is different from the lo and hi physical co-ordinates of the simulation domain.
    if (use_warpxba == false) dmap = amrex::DistributionMapping{ba};
    // Allocate output MultiFab for diagnostics. The data will be stored at cell-centers.
    int ngrow = (m_format == "sensei" || m_format == "ascent") ? 1 : 0;
    // The zero is hard-coded since the number of output buffers = 1 for FullDiagnostics
    m_mf_output[i_buffer][lev] = amrex::MultiFab(ba, dmap, m_varnames.size(), ngrow);


    if (lev == 0) {
        // The extent of the domain covered by the diag multifab, m_mf_output
        //default non-periodic geometry for diags
        amrex::Vector<int> diag_periodicity(AMREX_SPACEDIM,0);
        // Box covering the extent of the user-defined diagnostic domain
        amrex::Box domain = ba.minimalBox();
        // define geom object
        m_geom_output[i_buffer][lev].define( domain, &diag_dom,
                                             amrex::CoordSys::cartesian,
                                             diag_periodicity.data() );
    } else if (lev > 0) {
        // Take the geom object of previous level and refine it.
        m_geom_output[i_buffer][lev] = amrex::refine( m_geom_output[i_buffer][lev-1],
                                                      warpx.RefRatio(lev-1) );
    }
}


void
FullDiagnostics::InitializeFieldFunctors (int lev)
{
    auto & warpx = WarpX::GetInstance();

    // Clear any pre-existing vector to release stored data.
    m_all_field_functors[lev].clear();

    // Species index to loop over species that dump rho per species
    int i = 0;

    m_all_field_functors[lev].resize( m_varnames.size() );
    // Fill vector of functors for all components except individual cylindrical modes.
    for (int comp=0, n=m_all_field_functors[lev].size(); comp<n; comp++){
        if        ( m_varnames[comp] == "Ex" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 0), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Ey" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 1), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Ez" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 2), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Bx" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 0), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "By" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 1), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Bz" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 2), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "jx" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 0), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "jy" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 1), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "jz" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 2), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "rho" ){
            // Initialize rho functor to dump total rho
            m_all_field_functors[lev][comp] = std::make_unique<RhoFunctor>(lev, m_crse_ratio);
        } else if ( m_varnames[comp].rfind("rho_", 0) == 0 ){
            // Initialize rho functor to dump rho per species
            m_all_field_functors[lev][comp] = std::make_unique<RhoFunctor>(lev, m_crse_ratio, m_rho_per_species_index[i]);
            i++;
        } else if ( m_varnames[comp] == "F" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_F_fp(lev), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "G" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_G_fp(lev), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "phi" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_phi_fp(lev), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "part_per_cell" ){
            m_all_field_functors[lev][comp] = std::make_unique<PartPerCellFunctor>(nullptr, lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "part_per_grid" ){
            m_all_field_functors[lev][comp] = std::make_unique<PartPerGridFunctor>(nullptr, lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "divB" ){
            m_all_field_functors[lev][comp] = std::make_unique<DivBFunctor>(warpx.get_array_Bfield_aux(lev), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "divE" ){
            m_all_field_functors[lev][comp] = std::make_unique<DivEFunctor>(warpx.get_array_Efield_aux(lev), lev, m_crse_ratio);
        }
        else {
            amrex::Abort("Error: " + m_varnames[comp] + " is not a known field output type");
        }
    }
    AddRZModesToDiags( lev );
}


void
FullDiagnostics::PrepareFieldDataForOutput ()
{
    // First, make sure all guard cells are properly filled
    // Probably overkill/unnecessary, but safe and shouldn't happen often !!
    auto & warpx = WarpX::GetInstance();
    warpx.FillBoundaryE(warpx.getngE());
    warpx.FillBoundaryB(warpx.getngE());
    warpx.UpdateAuxilaryData();
    warpx.FillBoundaryAux(warpx.getngUpdateAux());

    // Update the RealBox used for the geometry filter in particle diags
    for (int i = 0; i < m_output_species.size(); ++i) {
        m_output_species[i].m_diag_domain = m_geom_output[0][0].ProbDomain();
    }
}

void
FullDiagnostics::MovingWindowAndGalileanDomainShift (int step)
{
    auto & warpx = WarpX::GetInstance();

    // Account for galilean shift
    amrex::Real new_lo[AMREX_SPACEDIM];
    amrex::Real new_hi[AMREX_SPACEDIM];
    const amrex::Real* current_lo = m_geom_output[0][0].ProbLo();
    const amrex::Real* current_hi = m_geom_output[0][0].ProbHi();

#if (AMREX_SPACEDIM == 3 )
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        new_lo[idim] = current_lo[idim] + warpx.m_galilean_shift[idim];
        new_hi[idim] = current_hi[idim] + warpx.m_galilean_shift[idim];
    }
#elif (AMREX_SPACEDIM == 2 )
    {
        new_lo[0] = current_lo[0] + warpx.m_galilean_shift[0];
        new_hi[0] = current_hi[0] + warpx.m_galilean_shift[0];
        new_lo[1] = current_lo[1] + warpx.m_galilean_shift[2];
        new_hi[1] = current_hi[1] + warpx.m_galilean_shift[2];
    }
#endif
    // Update RealBox of geometry with galilean-shifted boundary
    for (int lev = 0; lev < nmax_lev; ++lev) {
        m_geom_output[0][lev].ProbDomain( amrex::RealBox(new_lo, new_hi) );
    }

    // For Moving Window Shift
    if (warpx.moving_window_active(step+1)) {
        int moving_dir = warpx.moving_window_dir;
        amrex::Real moving_window_x = warpx.getmoving_window_x();
        // Get the updated lo and hi of the geom domain
        const amrex::Real* cur_lo = m_geom_output[0][0].ProbLo();
        const amrex::Real* cur_hi = m_geom_output[0][0].ProbHi();
        const amrex::Real* geom_dx = m_geom_output[0][0].CellSize();
        int num_shift_base = static_cast<int>((moving_window_x - cur_lo[moving_dir])
                                              / geom_dx[moving_dir]);
        // Update the diagnostic geom domain. Note that this is done only for the
        // base level 0 because m_geom_output[0][lev] share the same static RealBox
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            new_lo[idim] = cur_lo[idim];
            new_hi[idim] = cur_hi[idim];
        }
        new_lo[moving_dir] = cur_lo[moving_dir] + num_shift_base*geom_dx[moving_dir];
        new_hi[moving_dir] = cur_hi[moving_dir] + num_shift_base*geom_dx[moving_dir];
        // Update RealBox of geometry with shifted domain geometry for moving-window
        for (int lev = 0; lev < nmax_lev; ++lev) {
            m_geom_output[0][lev].ProbDomain( amrex::RealBox(new_lo, new_hi) );
        }
    }


}
