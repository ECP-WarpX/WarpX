#include "FullDiagnostics.H"
#include "WarpX.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/PartPerCellFunctor.H"
#include "ComputeDiagFunctors/PartPerGridFunctor.H"
#include "ComputeDiagFunctors/DivBFunctor.H"
#include "ComputeDiagFunctors/DivEFunctor.H"
#include "FlushFormats/FlushFormat.H"
#include "FlushFormats/FlushFormatPlotfile.H"
#include "FlushFormats/FlushFormatCheckpoint.H"
#include "FlushFormats/FlushFormatAscent.H"
#ifdef WARPX_USE_OPENPMD
#    include "FlushFormats/FlushFormatOpenPMD.H"
#endif
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
using namespace amrex::literals;


FullDiagnostics::FullDiagnostics (int i, std::string name)
    : Diagnostics(i, name)
{
    ReadParameters();
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
    if (m_species_names.size() == 0) m_species_names = mpc.GetSpeciesNames();
    // Initialize one ParticleDiag per species requested
    for (int i=0; i<m_species_names.size(); i++){
        const int idx = mpc.getSpeciesID(m_species_names[i]);
        m_all_species.push_back(ParticleDiag(m_diag_name, m_species_names[i],
                                             mpc.GetParticleContainerPtr(idx)));
    }
}

void
FullDiagnostics::ReadParameters ()
{
    // Read list of full diagnostics fields requested by the user.
    bool checkpoint_compatibility = BaseReadParameters();
    auto & warpx = WarpX::GetInstance();
    amrex::ParmParse pp(m_diag_name);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "plotfile" || m_format == "openpmd" ||
        m_format == "checkpoint" || m_format == "ascent" ||
        m_format == "sensei",
        "<diag>.format must be plotfile or openpmd or checkpoint or ascent or sensei");
    std::string period_string = "0";
    pp.query("period", period_string);
    m_intervals = IntervalsParser(period_string);
    bool raw_specified = pp.query("plot_raw_fields", m_plot_raw_fields);
    raw_specified += pp.query("plot_raw_fields_guards", m_plot_raw_fields_guards);

#ifdef WARPX_DIM_RZ
    pp.query("dump_rz_modes", m_dump_rz_modes);
#endif

    if (m_format == "checkpoint"){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            raw_specified == false &&
            checkpoint_compatibility == true,
            "For a checkpoint output, cannot specify these parameters as all data must be dumped "
            "to file for a restart");
    }
    // Number of buffers = 1 for FullDiagnostics.
    // It is used to allocate the number of output multu-level MultiFab, m_mf_output
    m_num_buffers = 1;
}

void
FullDiagnostics::Flush ( int i_buffer )
{
    // This function should be moved to Diagnostics when plotfiles/openpmd format
    // is supported for BackTransformed Diagnostics, in BTDiagnostics class.
    auto & warpx = WarpX::GetInstance();
    m_flush_format->WriteToFile(
        m_varnames, m_mf_output[i_buffer], warpx.Geom(), warpx.getistep(),
        warpx.gett_new(0), m_all_species, nlev_output, m_file_prefix,
        m_plot_raw_fields, m_plot_raw_fields_guards, m_plot_raw_rho, m_plot_raw_F);

    FlushRaw();
}

void
FullDiagnostics::FlushRaw () {}


bool
FullDiagnostics::DoDump (int step, int i_buffer, bool force_flush)
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

    // First index of m_all_field_functors[lev] where RZ modes are stored
    int icomp = m_all_field_functors[0].size();
    const std::array<std::string, 3> coord {"r", "theta", "z"};

    // Er, Etheta, Ez, Br, Btheta, Bz, jr, jtheta, jz
    // Each of them being a multi-component multifab
    int n_new_fields = 9;
    if (divE_requested) {
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
    // Sum the number of components in input vector m_all_field_functors
    // and check that it corresponds to the number of components in m_varnames
    // and m_mf_output
    int ncomp_from_src = 0;
    for (int i=0; i<m_all_field_functors[0].size(); i++){
        ncomp_from_src += m_all_field_functors[lev][i]->nComp();
    }
    AMREX_ALWAYS_ASSERT( ncomp_from_src == m_varnames.size() );
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE (
        ba.coarsenable(m_crse_ratio),
        "Invalid coarsening ratio for warpx boxArray. Must be an integer divisor of the blocking factor."
    );

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
        // blocking factor at level, lev.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( blockingFactor[idim] % m_crse_ratio[idim]==0,
                       " coarsening ratio must be integer divisor of blocking factor");
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
            m_lo[idim] = warpx.Geom(lev).ProbLo(idim) + warpx.Geom(lev).CellSize(idim)/2.0_rt +
                ba.getCellCenteredBox(0).smallEnd(idim) * warpx.Geom(lev).CellSize(idim);
            m_hi[idim] = warpx.Geom(lev).ProbLo(idim) + warpx.Geom(lev).CellSize(idim)/2.0_rt +
                ba.getCellCenteredBox( ba.size()-1 ).bigEnd(idim) * warpx.Geom(lev).CellSize(idim);
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
    int ngrow = (m_format == "sensei") ? 1 : 0;
    // The zero is hard-coded since the number of output buffers = 1 for FullDiagnostics
    m_mf_output[i_buffer][lev] = amrex::MultiFab(ba, dmap, m_varnames.size(), ngrow);

}


void
FullDiagnostics::InitializeFieldFunctors (int lev)
{
    auto & warpx = WarpX::GetInstance();
    // Clear any pre-existing vector to release stored data.
    m_all_field_functors[lev].clear();

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
            // rho_new is stored in component 1 of rho_fp when using PSATD
#ifdef WARPX_USE_PSATD
            amrex::MultiFab* rho_new = new amrex::MultiFab(*warpx.get_pointer_rho_fp(lev), amrex::make_alias, 1, 1);
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(rho_new, lev, m_crse_ratio);
#else
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_rho_fp(lev), lev, m_crse_ratio);
#endif
        } else if ( m_varnames[comp] == "F" ){
            m_all_field_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_F_fp(lev), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "part_per_cell" ){
            m_all_field_functors[lev][comp] = std::make_unique<PartPerCellFunctor>(nullptr, lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "part_per_grid" ){
            m_all_field_functors[lev][comp] = std::make_unique<PartPerGridFunctor>(nullptr, lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "divB" ){
            m_all_field_functors[lev][comp] = std::make_unique<DivBFunctor>(warpx.get_array_Bfield_aux(lev), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "divE" ){
            m_all_field_functors[lev][comp] = std::make_unique<DivEFunctor>(warpx.get_array_Efield_aux(lev), lev, m_crse_ratio);
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
    warpx.FillBoundaryE(warpx.getngE(), warpx.getngExtra());
    warpx.FillBoundaryB(warpx.getngE(), warpx.getngExtra());
#ifndef WARPX_USE_PSATD
    warpx.FillBoundaryAux(warpx.getngUpdateAux());
#endif
    warpx.UpdateAuxilaryData();
}
