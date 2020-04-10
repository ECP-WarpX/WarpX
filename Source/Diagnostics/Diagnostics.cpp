#include "Diagnostics.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/PartPerCellFunctor.H"
#include "ComputeDiagFunctors/PartPerGridFunctor.H"
#include "ComputeDiagFunctors/DivBFunctor.H"
#include "ComputeDiagFunctors/DivEFunctor.H"
#include "WarpX.H"
#include "Utils/Average.H"
#include "Utils/WarpXUtil.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : diag_name(name), diag_index(i)
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
    file_prefix = "diags/" + diag_name;
    pp.query("file_prefix", file_prefix);
    pp.query("period", m_period);
    pp.query("plot_raw_fields", m_plot_raw_fields);
    pp.query("plot_raw_fields_guards", m_plot_raw_fields_guards);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_plot_F==false, "cannot plot_F yet");
    if (!pp.queryarr("fields_to_plot", varnames)){
        varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    }
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


    //// Read user-defined (lo,hi) extents for the diag
    m_lo.resize(AMREX_SPACEDIM);
    m_hi.resize(AMREX_SPACEDIM);

    if (!pp.queryarr("diag_lo", m_lo, 0, AMREX_SPACEDIM) ) {
       for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
            m_lo[idim] = warpx.Geom(0).ProbLo(idim);
       }
    }
    if (! pp.queryarr("diag_hi", m_hi, 0, AMREX_SPACEDIM) ) {
       for (int idim =0; idim < AMREX_SPACEDIM; ++idim) {
            m_hi[idim] = warpx.Geom(0).ProbHi(idim);
       }
    }

    // Initialize cr_ratio with default value of 1 for each dimension.
    Vector<int> cr_ratio(AMREX_SPACEDIM, 1);
    // Read user-defined coarsening ratio for the output multifab.
    if (pp.queryarr("coarsening_ratio", cr_ratio, 0, AMREX_SPACEDIM) ) {
       for (int idim =0; idim < AMREX_SPACEDIM; ++idim) {
           m_crse_ratio[idim] = cr_ratio[idim];
       }
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
        all_field_functors[lev].resize( varnames.size() );
        // Fill vector of functors for all components except individual
        // cylindrical modes
        for (int comp=0, n=all_field_functors[lev].size(); comp<n; comp++){
            if        ( varnames[comp] == "Ex" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 0), lev, m_crse_ratio);
            } else if ( varnames[comp] == "Ey" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 1), lev, m_crse_ratio);
            } else if ( varnames[comp] == "Ez" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, 2), lev, m_crse_ratio);
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 0), lev, m_crse_ratio);
            } else if ( varnames[comp] == "By" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 1), lev, m_crse_ratio);
            } else if ( varnames[comp] == "Bz" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, 2), lev, m_crse_ratio);
            } else if ( varnames[comp] == "jx" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 0), lev, m_crse_ratio);
            } else if ( varnames[comp] == "jy" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 1), lev, m_crse_ratio);
            } else if ( varnames[comp] == "jz" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_current_fp(lev, 2), lev, m_crse_ratio);
            } else if ( varnames[comp] == "rho" ){
                // rho_new is stored in component 1 of rho_fp when using PSATD
#ifdef WARPX_USE_PSATD
                MultiFab* rho_new = new MultiFab(*warpx.get_pointer_rho_fp(lev), amrex::make_alias, 1, 1);
                all_field_functors[lev][comp] = new CellCenterFunctor(rho_new, lev, m_crse_ratio);
#else
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_rho_fp(lev), lev, m_crse_ratio);
#endif
            } else if ( varnames[comp] == "F" ){
                all_field_functors[lev][comp] = new CellCenterFunctor(warpx.get_pointer_F_fp(lev), lev, m_crse_ratio);
            } else if ( varnames[comp] == "part_per_cell" ){
                all_field_functors[lev][comp] = new PartPerCellFunctor(nullptr, lev, m_crse_ratio);
            } else if ( varnames[comp] == "part_per_grid" ){
                all_field_functors[lev][comp] = new PartPerGridFunctor(nullptr, lev, m_crse_ratio);
            } else if ( varnames[comp] == "divB" ){
                all_field_functors[lev][comp] = new DivBFunctor(warpx.get_array_Bfield_aux(lev), lev, m_crse_ratio);
            } else if ( varnames[comp] == "divE" ){
                all_field_functors[lev][comp] = new DivEFunctor(warpx.get_array_Efield_aux(lev), lev, m_crse_ratio);
            }
        }

        AddRZModesToDiags( lev );

        // At this point, varnames.size() >= all_field_functors[0].size()
        DefineDiagMultiFab( lev );

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
    int icomp_dst = 0;
    for(int lev=0; lev<nlev; lev++){
        for (int icomp=0, n=all_field_functors[0].size(); icomp<n; icomp++){
            // Call all functors in all_field_functors[lev]. Each of them computes
            // a diagnostics and writes in one or more components of the output
            // multifab mf_avg[lev].
            all_field_functors[lev][icomp]->operator()(mf_avg[lev], icomp_dst);
            // update the index of the next component to fill
            icomp_dst += all_field_functors[lev][icomp]->nComp();
        }
    }
    // Check that the proper number of components of mf_avg were updated.
    AMREX_ALWAYS_ASSERT( icomp_dst == varnames.size() );
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
Diagnostics::AddRZModesToDiags (int lev)
{
#ifdef WARPX_DIM_RZ
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

    // First index of all_field_functors[lev] where RZ modes are stored
    int icomp = all_field_functors[0].size();
    const std::array<std::string, 3> coord {"r", "theta", "z"};

    // Er, Etheta, Ez, Br, Btheta, Bz, jr, jtheta, jz
    // Each of them being a multi-component multifab
    all_field_functors[lev].resize( all_field_functors[0].size() + 9 );
    // E
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        all_field_functors[lev][icomp] = new
            CellCenterFunctor(warpx.get_pointer_Efield_aux(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        AddRZModesToOutputNames(std::string("E") + coord[dim],
                                warpx.get_pointer_Efield_aux(0, 0)->nComp());
        icomp += 1;
    }
    // B
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        all_field_functors[lev][icomp] = new
            CellCenterFunctor(warpx.get_pointer_Bfield_aux(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        AddRZModesToOutputNames(std::string("B") + coord[dim],
                                warpx.get_pointer_Bfield_aux(0, 0)->nComp());
        icomp += 1;
    }
    // j
    for (int dim=0; dim<3; dim++){
        // 3 components, r theta z
        all_field_functors[lev][icomp] = new
            CellCenterFunctor(warpx.get_pointer_current_fp(lev, dim), lev,
                              m_crse_ratio, false, ncomp_multimodefab);
        icomp += 1;
        AddRZModesToOutputNames(std::string("J") + coord[dim],
                                warpx.get_pointer_current_fp(0, 0)->nComp());
    }
    // Sum the number of components in input vector all_field_functors
    // and check that it corresponds to the number of components in varnames
    // and mf_avg
    int ncomp_from_src = 0;
    for (int i=0; i<all_field_functors[0].size(); i++){
        ncomp_from_src += all_field_functors[lev][i]->nComp();
    }
    AMREX_ALWAYS_ASSERT( ncomp_from_src == varnames.size() );
#endif
}

void
Diagnostics::AddRZModesToOutputNames (const std::string& field, int ncomp){
#ifdef WARPX_DIM_RZ
    // In cylindrical geometry, real and imag part of each mode are also
    // dumped to file separately, so they need to be added to varnames
    varnames.push_back( field + "_0_real" );
    for (int ic=1 ; ic < ncomp ; ic += 2) {
        varnames.push_back( field + "_" + std::to_string(ic) + "_real" );
        varnames.push_back( field + "_" + std::to_string(ic) + "_imag" );
    }
#endif
}

void
Diagnostics::DefineDiagMultiFab ( int lev ) {

    auto & warpx = WarpX::GetInstance();
    amrex::RealBox diag_dom;
    bool use_warpxba = true;

    // Default BoxArray and DistributionMap for mf_avg from the warpx instance
    BoxArray ba = warpx.boxArray(lev);
    DistributionMapping dmap = warpx.DistributionMap(lev);

    // Find if user-defined diag lo and hi are different from the simulation
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
         // To ensure that the diag lo and hi are within the domain defined at lev
         diag_dom.setLo(idim, max(m_lo[idim],warpx.Geom(lev).ProbLo(idim)) );
         diag_dom.setHi(idim, min(m_hi[idim],warpx.Geom(lev).ProbHi(idim)) );
         if ( fabs(warpx.Geom(lev).ProbLo(idim) - diag_dom.lo(idim))
                                >  warpx.Geom(lev).CellSize(idim) )
              use_warpxba = false;
         if ( fabs(warpx.Geom(lev).ProbHi(idim) - diag_dom.hi(idim))
                                > warpx.Geom(lev).CellSize(idim) )
              use_warpxba = false;
    }

    if (use_warpxba == false) {
       // Following are the steps to create a new box array
       IntVect lo(AMREX_D_DECL(0,0,0));
       IntVect hi(AMREX_D_DECL(1,1,1));
       for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
           // lo index
           lo[idim] = static_cast<int>( floor (
                                ( diag_dom.lo(idim)
                                - warpx.Geom(lev).ProbLo(idim) )
                                / warpx.Geom(lev).CellSize(idim) ) );
           // hi index
           hi[idim] = static_cast<int> ( ceil (
                                ( diag_dom.hi(idim)
                                - warpx.Geom(lev).ProbLo(idim) )
                                / warpx.Geom(lev).CellSize(idim) ) ) ;

           // Modify lo and hi if the crse ratio is power of two
           // But the lo and hi are still not coarsenable
           int mod_lo = lo[idim] % m_crse_ratio[idim];
           int mod_hi = hi[idim] % m_crse_ratio[idim];
           if (mod_lo > 0) lo[idim] -= mod_lo;
           if (mod_hi > 0) hi[idim] += (m_crse_ratio[idim] - mod_hi);
           // Subtract hi_index by 1 accounting for cell-centered indextype of mf_avg
           hi[idim]  -= 1;
           // if hi<=lo, then hi = lo + 1, to ensure one cell in that dimension
           if ( hi[idim] <= lo[idim] ) {
                hi[idim]  = lo[idim] + 1;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                   m_crse_ratio[idim]==1, "coarsening ratio in reduced dimension must be 1."
                );
           }
           // update lo and hi according to the modified ilo and ihi in index space
           m_lo[idim] = lo[idim] * warpx.Geom(lev).CellSize(idim)
                      + warpx.Geom(lev).ProbLo(idim);
           m_hi[idim] = hi[idim] * warpx.Geom(lev).CellSize(idim)
                      + warpx.Geom(lev).ProbLo(idim);
       }
       // Box for the reduced-domain diag
       Box diag_box( lo, hi );
       // Define box array
       BoxArray diag_ba;
       diag_ba.define(diag_box);
       ba = diag_ba.maxSize( warpx.maxGridSize( lev ) );
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_crse_ratio.min() > 0,
        "Coarsening ratio must be non-zero."
    );
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE (
        ba.coarsenable(m_crse_ratio),
        "Invalid coarsening ratio for diagnostics. Must be a power of two and an integer divisor of the blocking factor."
    );
    // The boxArray is coarsened based on the user-defined coarsening ratio
    ba.coarsen(m_crse_ratio);
    // Generate a new distribution map if lo and hi for diag is different from simulation
    if (use_warpxba == false) dmap = DistributionMapping{ba};
    // Allocate output multifab
    // Note: default MultiFab constructor is cell-centered
    mf_avg[lev] = MultiFab(ba, dmap, varnames.size(), 0);

}
